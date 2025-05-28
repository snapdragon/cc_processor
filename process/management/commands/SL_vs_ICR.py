import logging
from time import sleep
import json
from gprofiler import GProfiler
import pandas as pd
import requests
import io
import csv
from collections import defaultdict

from django.core.management.base import BaseCommand

from process.models import Project, Protein, Phospho, Abundance, Statistic

from process.constants import ABUNDANCES_RAW, ABUNDANCES_NORMALISED_LOG2_MEAN, ANOVA, Q_VALUE, CURVE_FOLD_CHANGE, PROTEIN_MAX_Q

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    def handle(self, *args, **options):
        logger.info("Compare results for SL and ICR")

        SL_project = Project.objects.get(name="SL")
        ICR_project = Project.objects.get(name="ICR")
        Original_project = Project.objects.get(name="Original")

        # print_percentages(project)

        print_thresholds(SL_project)
        print_thresholds(ICR_project)
        print_thresholds(Original_project)



def print_percentages(project):
    # Output percentages of blank cells for proteins and phosphos

    protein_rows = Protein.objects.filter(project=project).count()
    phospho_rows = Phospho.objects.filter(protein__project=project).count()

    raw_protein_abundances_num = Abundance.objects.filter(
        statistic__protein__project = project,
        statistic__statistic_type__name = ABUNDANCES_RAW,
        replicate__mean = False
    ).count()

    raw_phospho_abundances_num = Abundance.objects.filter(
        statistic__phospho__protein__project = project,
        statistic__statistic_type__name = ABUNDANCES_RAW,
        replicate__mean = False
    ).count()

    protein_blank_percentage = round((raw_protein_abundances_num / (protein_rows * 30)) * 100, 2)
    phospho_blank_percentage = round((raw_phospho_abundances_num / (phospho_rows * 30)) * 100, 2)

    print(f"Total populated proteins: {raw_protein_abundances_num} of {protein_rows * 30}: {protein_blank_percentage}% populated")
    print(f"Total populated phosphos: {raw_phospho_abundances_num} of {phospho_rows * 30}: {phospho_blank_percentage}% populated")




def print_thresholds(project):
    # Output percentages of proteins above and below q_value and curve_fold_change limits

    protein_rows = Protein.objects.filter(project=project).count()

    statistics = Statistic.objects.filter(
        statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
        protein__project=project,
    )

    num_under_q_value = 0
    num_over_curve_fold = 0

    for statistic in statistics:
        if (
            statistic.metrics is None
            or statistic.metrics.get(ANOVA) is None
            or statistic.metrics[ANOVA].get(Q_VALUE) is None
            or statistic.metrics[ANOVA][Q_VALUE] >= 0.01
        ):
            continue

        num_under_q_value += 1

    for statistic in statistics:
        if (
            statistic.metrics.get(CURVE_FOLD_CHANGE) is None
            or round(statistic.metrics[CURVE_FOLD_CHANGE], 1) < 1.2
        ):
            continue

        num_over_curve_fold += 1

    q_value_percentage = round((num_under_q_value / len(statistics)) * 100, 2)
    curve_fold_percentage = round((num_over_curve_fold / len(statistics)) * 100, 2)

    print(f"Project: {project.name}")

    print(f"Total under q value: {num_under_q_value} of {len(statistics)}: {q_value_percentage}%")
    print(f"Total over curve fold: {num_over_curve_fold} of {len(statistics)}: {curve_fold_percentage}%")

    print("\n")
