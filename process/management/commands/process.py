import logging
import statistics

from django.core.management.base import BaseCommand
from django.db.models.query import QuerySet

from process.models import ColumnName, Project, ProteinReading, Replicate

# TODO - make this configurable by flag
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

PROJECT_NAME = "Soliman Labs"


class Command(BaseCommand):
    help = "Processes all proteins for a given project"

    def handle(self, *args, **kwargs):
        # TODO - make this a command line option, to cater for other projects
        project = Project.objects.get(name=PROJECT_NAME)
        replicates = Replicate.objects.filter(project=project)
        protein_readings = ProteinReading.objects.filter(
            column_name__replicate__project=project
        )

        self._process(project, replicates, protein_readings)

    def _process(self, project, replicates, protein_readings):
        """
        Does all the required calculations. The steps are:

        TODO - put the descriptions below in the function comments, not here, along with
            data structures.

        1) Calculate the median for each replicate for each stage.
            In other words, get all the values for each column and find the median.

        2) Calculate the mean for each protein for each stage across replicates.
            So for proteins 'ABCD', 'EFGH' with replicates 'One', 'Two' and stages '1h', '2h',
            it will get the abundances for 'ABCD 1h' for each of the replicates, then take
            the mean.

        3) TODO
        """

        # TODO - make each call flaggable
        medians = self._all_replicates(
            func=self._calc_replicate_column_medians,
            replicates=replicates,
            protein_readings=protein_readings,
        )

        print(medians)

        abundance_means_across_replicates_by_stage = (
            self._calculate_abundance_means_across_replicates_by_stage(protein_readings)
        )

        print(abundance_means_across_replicates_by_stage)

    def _calculate_abundance_means_across_replicates_by_stage(
        self, protein_readings: QuerySet[ProteinReading]
    ):
        logger.info("Calculating mean across replicates by stage for each protein")

        means: dict = {}

        protein_no = 0

        for protein_reading in protein_readings:
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Appending for {protein_no}, {protein_reading.protein.accession_number} {protein_reading.column_name.sample_stage.name}",
            )

            protein = protein_reading.protein

            if not means.get(protein):
                means[protein] = {}

            stage = protein_reading.column_name.sample_stage

            if not means[protein].get(stage):
                means[protein][stage] = []

            if protein_reading.reading:
                means[protein][stage].append(protein_reading.reading)

        protein_stage_no = 0

        for protein in means.keys():
            for stage in means[protein]:
                protein_stage_no += 1

                abundances = means[protein][stage]

                if len(abundances):
                    mean = sum(abundances) / len(abundances)
                    means[protein][stage] = mean
                else:
                    means[protein][stage] = None

                protein_stage_no += 1
                self._count_logger(
                    protein_stage_no,
                    10000,
                    f"Mean for {protein_stage_no}, {protein} {stage}: {mean}",
                )

        return means

    def _all_replicates(self, *args, **kwargs):
        """
        Calls the passed function for each replicate for the project.
        """
        results = {}

        # Remove the passed function and replicates as they're not needed by the passed function
        call_kwargs = kwargs.copy()
        func = call_kwargs.pop("func")
        replicates = call_kwargs.pop("replicates")

        for replicate in replicates:
            call_kwargs["replicate"] = replicate

            results[replicate] = func(**call_kwargs)

        return results

    def _calc_replicate_column_medians(
        self, replicate: Replicate, protein_readings: QuerySet[ProteinReading]
    ):
        logger.info("Calculating column medians by replicate")

        column_medians = {}

        column_names = ColumnName.objects.filter(replicate__name=replicate.name)

        for column_name in column_names:
            readings = []

            protein_readings_by_column = protein_readings.filter(
                column_name=column_name
            )

            for protein_reading in protein_readings_by_column:
                if protein_reading.reading:
                    readings.append(protein_reading.reading)
                else:
                    # TODO - what to do about None values?
                    readings.append(0)

            median = statistics.median(readings)

            logger.info(
                f"calculating median for replicate {replicate.name} stage {column_name.sample_stage.name}: {median}"
            )

            column_medians[column_name.id] = median

        return column_medians

    def _count_logger(self, i: int, step: int, output: str):
        if i % step == 0:
            logger.info(output)
