import statistics

from django.core.management.base import BaseCommand

from process.models import ColumnName, ProteinReading

REPLICATE_NAMES = ("One", "Two", "Three")


class Command(BaseCommand):
    help = "Runs the custom script"

    def handle(self, *args, **kwargs):
        self._calc_replicate_column_medians(REPLICATE_NAMES[0])
        self._calc_replicate_column_medians(REPLICATE_NAMES[1])
        self._calc_replicate_column_medians(REPLICATE_NAMES[2])

    def _calc_replicate_column_medians(self, replicate_name: str):
        column_medians = {}

        column_names = ColumnName.objects.filter(replicate__name=replicate_name)

        for column_name in column_names:
            readings = []

            protein_readings = ProteinReading.objects.filter(column_name=column_name)

            for protein_reading in protein_readings:
                # NaN does not equal NaN
                if protein_reading.reading != protein_reading.reading:
                    # TODO - what to do about NaN values?
                    readings.append(0)
                else:
                    readings.append(protein_reading.reading)

            median = statistics.median(readings)

            column_medians[column_name.id] = median

        return column_medians
