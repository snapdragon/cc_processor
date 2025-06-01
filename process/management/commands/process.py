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
        # TODO - make each call flaggable
        medians = self._all_replicates(
            func=self._calc_replicate_column_medians,
            replicates=replicates,
            protein_readings=protein_readings,
        )

        print(medians)

        # protein_abundance_means_by_stage = self._calculate_abundance_means_across_replicates(protein_readings)

    # def _calculate_abundance_means_across_replicates(self, protein_readings: QuerySet[ProteinReading]):
    #     means = {}

    #     for protein_reading in protein_readings:
    #         protein = protein_reading.protein

    #         if not means.get(protein):
    #             means[protein] = {}

    #         stage = protein_reading.column_name.sample_stage

    #         if not means[protein].get(stage):
    #             means[protein][stage] = []

    #         means[protein][stage].append(protein_reading.reading)

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
        column_medians = {}

        column_names = ColumnName.objects.filter(replicate__name=replicate.name)

        for column_name in column_names:
            readings = []

            protein_readings_by_column = protein_readings.filter(
                column_name=column_name
            )

            for protein_reading in protein_readings_by_column:
                # NaN does not equal NaN
                if protein_reading.reading != protein_reading.reading:
                    # TODO - what to do about NaN values?
                    readings.append(0)
                else:
                    readings.append(protein_reading.reading)

            median = statistics.median(readings)

            logger.info(
                f"calculating median for replicate {replicate.name} stage {column_name.sample_stage.name}: {median}"
            )

            column_medians[column_name.id] = median

        return column_medians
