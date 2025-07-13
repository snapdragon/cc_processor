# Look out for 'No Original ANOVA' errors

import logging
import json
import math

import pandas as pd
from django.core.management.base import BaseCommand

from process.models import (
    Project,
    Protein,
    Abundance,
    Phospho,
    StatisticType,
    Statistic,
    Replicate,
    SampleStage
)

from process.constants import (
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    P_VALUE,
    F_STATISTICS,
    Q_VALUE,
    ANOVA,
    G_STATISTIC,
    FISHER_G,
    FREQUENCY,
    CURVE_FOLD_CHANGE,
)

# TODO - move to constants file
METRICS_STRING_FIELDS = ['peak_average', 'curve_peak']

METRICS_NUMBER_FIELDS = [
    'standard_deviation',
    'variance_average',
    'skewness_average',
    'kurtosis_average',
    'max_fold_change_average',
    'residuals_average',
    'R_squared_average',
    'residuals_all',
    'R_squared_all',
    'curve_fold_change'
]

METRICS_ANOVA_FIELDS = [
    P_VALUE,
    F_STATISTICS,
    Q_VALUE,
]

METRICS_FISHER_G_FIELDS = [
    G_STATISTIC,
    P_VALUE,
    FREQUENCY,
    Q_VALUE,
]

TOLERANCE = 0.1

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to import",
        )
        parser.add_argument(
            "--accession-number",
            required=False,
            help="The accession number of the protein to compare.",
        )

    def handle(self, *args, **options):
        logger.info("Listing proteins and phosphos of interest")

        accession_number = options["accession_number"]
        project_name = options["project"]

        if not project_name:
            raise Exception(f"Invalud project name {project_name}")

        project = Project.objects.get(name=project_name)

        if accession_number:
            statistics = Statistic.objects.filter(
                statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                protein__project=project,
                protein__accession_number=accession_number
            )
        else:
            statistics = Statistic.objects.filter(
                statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                protein__project=project
            )

        logger.info("Proteins of interest")

        self.run_stats(statistics, 0.02)

        if accession_number:
            statistics = Statistic.objects.filter(
                statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                phospho__protein__project=project,
                phospho__protein__accession_number=accession_number
            )
        else:
            statistics = Statistic.objects.filter(
                statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                phospho__protein__project=project
            )

        logger.info("Phosphos of interest")

        self.run_stats(statistics, 0.02)


    def run_stats(self, statistics, max_q):
        total = 0
        num_matches = 0

        for statistic in statistics:
            total += 1

            match = True

            if statistic.metrics[ANOVA][Q_VALUE] >= max_q:
                match = False

            if statistic.metrics.get(CURVE_FOLD_CHANGE) is None or statistic.metrics[CURVE_FOLD_CHANGE] <= 1.2:
                match = False

            if match:
                num_matches += 1
                print(f"Match for {statistic.phospho.protein.accession_number}")

        print(f"Total protein matches {num_matches} of {total}")
