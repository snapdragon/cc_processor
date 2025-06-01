import statistics

from django.core.management.base import BaseCommand
from django.db.models.query import QuerySet

from process.models import ColumnName, Project, ProteinReading, Replicate

PROJECT_NAME = "Soliman Labs"


class Command(BaseCommand):
    help = "Runs the custom script"

    def handle(self, *args, **kwargs):
        # TODO - make this a command line option, to cater for other projects
        project = Project.objects.get(name=PROJECT_NAME)

        replicates = Replicate.objects.filter(project=project)
        protein_readings = ProteinReading.objects.filter(
            column_name__replicate__project=project
        )

        medians = self._all_replicates(
            self._calc_replicate_column_medians, replicates, protein_readings
        )

        print(medians)

    def _all_replicates(self, func, replicates, protein_readings):
        returns = {}

        for replicate in replicates:
            returns[replicate] = func(replicate, protein_readings)

        return returns

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

            print(
                f"calculating median for replicate {replicate.name} stage {column_name.sample_stage.name}: {median}"
            )

            column_medians[column_name.id] = median

        return column_medians
