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
    ABUNDANCES_RAW,
    RAW,
    PROTEIN_ABUNDANCES,
    ABUNDANCES_NORMALISED_MEDIAN,
    ABUNDANCES_NORMALISED_LOG2_ARREST,
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ABUNDANCES_NORMALISED_MIN_MAX,
    ABUNDANCES_NORMALISED_ZERO_MAX,
    ABUNDANCES_IMPUTED,
    P_VALUE,
    F_STATISTICS,
    Q_VALUE,
    ANOVA,
    G_STATISTIC,
    FISHER_G,
    FREQUENCY,
    PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
    PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX,
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
            "--accession-number",
            required=False,
            help="The accession number of the protein to compare.",
        )

    def handle(self, *args, **options):
        accession_number = options["accession_number"]

        logger.info(f"Comparing original and new")

        project = Project.objects.get(name="ICR")

        proteins = Protein.objects.filter(project=project)

        statistic_types = StatisticType.ojects.all()

        for statistic_type in statistic_types:
            num_abundances = Abundance.objects.filter(
                statistic__statistic_type=statistic_type,
            )

            print(f"{statistic_type.name} {num_abundances}")

