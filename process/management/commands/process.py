import logging
import statistics
from typing import Final

from django.core.management.base import BaseCommand
from django.db.models.query import QuerySet

from process.models import ColumnName, Project, ProteinReading, Replicate

# TODO - make this configurable by flag
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "Processes all proteins for a given project"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to process",
        )

    def handle(self, *args, **options):
        project_name = options["project"]

        logger.info("Processing for project {project_name}")

        project: Final = Project.objects.get(name=project_name)
        replicates: Final = Replicate.objects.filter(project=project)
        column_names: Final = ColumnName.objects.filter(replicate__project=project)
        protein_readings: Final = ProteinReading.objects.filter(
            column_name__replicate__project=project
        )

        self._process(project, replicates, protein_readings, column_names)

    def _process(self, project, replicates, protein_readings, column_names):
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

        3) Normalise all abundances. This is done by dividing each abundance by the median
            for its column, then averaging them across replicates.
        """

        # TODO - make each call flaggable
        # TODO - rename 'medians' to something more informative
        # TODO - does this need to be by replicate? Why not just all columns at once?
        medians = self._all_replicates(
            func=self._calc_replicate_column_medians,
            replicates=replicates,
            protein_readings=protein_readings,
            column_names=column_names,
        )

        print(f"Medians for project {project.name}")
        for replicate in medians.keys():
            print(f"Replicate: {replicate.name}")

            for stage in medians[replicate].keys():
                print(f"    {stage.name}: {medians[replicate][stage]}")

        return

        means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(protein_readings)
        )

        logger.info(means_across_replicates_by_stage)

        normalised_protein_readings = self._calculate_first_level_normalisation(
            protein_readings, medians
        )

        logger.info(
            f"Number of normalised readings: f{len(normalised_protein_readings)}"
        )

        normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                normalised_protein_readings
            )
        )

        logger.info(normalised_means_across_replicates_by_stage)

        arrest_log2_normalised_protein_readings = (
            self._calculate_arrest_log2_normalisation(protein_readings, medians)
        )

        logger.info(
            f"Number of arrest log2 normalised readings: f{len(normalised_protein_readings)}"
        )

        arrest_log2_normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                arrest_log2_normalised_protein_readings
            )
        )

        logger.info(arrest_log2_normalised_means_across_replicates_by_stage)

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

    def _calculate_arrest_log2_normalisation(
        self, normalised_protein_readings: QuerySet[ProteinReading]
    ):
        pass
        # logger.info("log2 arrest normalising abundances")

        # # TODO - this will need to be changed for per-project. Put it in the Project model maybe?
        # ARRESTING_AGENT = "Nocodozole"

        # protein_arrest_values = {}

        # # Get all Nocodozole values by replicate by protein
        # # TODO - this is super inefficient
        # for protein_reading in normalised_protein_readings:
        #     if protein_reading.column_name.sample_stage.name == ARRESTING_AGENT:
        #         protein_arrest_values[protein_reading.column_name.replicate] = {
        #             [protein_reading.protein] = protein_reading.reading

        # arrest_logs_normalised_protein_readings = []

        # for protein_reading in protein_readings:
        #     reading = protein_reading.reading
        #     normalised_reading = None

        #     if reading:
        #         # TODO - need to round this?
        #         normalised_reading = reading / protein_arrest_values[protein_reading.protein]

        #     # TODO - not a QuerySet
        #     # TODO - inefficient
        #     arrest_logs_normalised_protein_readings.append(ProteinReading(protein=protein_reading.protein, column_name = protein_reading.column_name, reading=normalised_reading))

        # return arrest_logs_normalised_protein_readings

    def _calculate_first_level_normalisation(
        self, protein_readings: QuerySet[ProteinReading], medians
    ):
        logger.info("Normalising abundances")

        normalised_protein_readings = []

        protein_no = 0

        for protein_reading in protein_readings:
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Normalising for {protein_no}, {protein_reading.protein.accession_number} {protein_reading.column_name.sample_stage.name}",
            )

            reading = protein_reading.reading

            normalised_reading = None

            if reading:
                median = medians[protein_reading.column_name.replicate][
                    protein_reading.column_name
                ]

                # TODO - need to round this?
                # TODO - might there be undefined medians? If so what to do?
                #   Maybe throw an exception in the median calculation func if one is encountered?
                #   They probably shouldn't happen.
                normalised_reading = reading / median

            # TODO - not a QuerySet
            # TODO - inefficient
            normalised_protein_readings.append(
                ProteinReading(
                    protein=protein_reading.protein,
                    column_name=protein_reading.column_name,
                    reading=normalised_reading,
                )
            )

        return normalised_protein_readings

    def _calculate_means_across_replicates_by_stage(
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
                f"Normalising for {protein_no}, {protein_reading.protein.accession_number} {protein_reading.column_name.sample_stage.name}",
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

    def _calc_replicate_column_medians(
        self,
        replicate: Replicate,
        protein_readings: QuerySet[ProteinReading],
        column_names: QuerySet[ColumnName],
    ):
        logger.info("Calculating column medians by replicate")

        column_medians = {}

        column_names_by_replicate = column_names.filter(replicate__name=replicate.name)

        for column_name in column_names_by_replicate:
            readings = []

            protein_readings_by_column = protein_readings.filter(
                column_name=column_name
            )

            for protein_reading in protein_readings_by_column:
                # TODO - what to do about None values?
                if protein_reading.reading is not None:
                    readings.append(protein_reading.reading)

            median = statistics.median(readings)

            column_medians[column_name] = median

        return column_medians

    def _count_logger(self, i: int, step: int, output: str):
        if i % step == 0:
            logger.info(output)
