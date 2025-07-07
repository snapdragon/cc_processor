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
    PROTEIN_ABUNDANCES_RAW,
    RAW,
    PROTEIN_ABUNDANCES,
    PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
    PROTEIN_ABUNDANCES_NORMALISED_LOG2_ARREST,
    PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
    PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX,
    PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX,
    PROTEIN_ABUNDANCES_IMPUTED,
    P_VALUE,
    F_STATISTICS,
    Q_VALUE,
    ANOVA,
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
    # Q_VALUE,
]

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

        project_original = Project.objects.get(name="Original")
        project_process = Project.objects.get(name="ICR")

        if not project_process.with_bugs:
            raise Exception("Only a with_bugs ICR project can be compared to the original. You'll have to run 'process' again with --with-bugs.")

        if accession_number:
            self._compare(accession_number, project_original, project_process)
        else:
            # TODO - only compares proteins in original, any missing from
            #   new will be missed. Do both?
            proteins = Protein.objects.filter(project=project_original)

            num_proteins = 0

            for pr in proteins:
                if not num_proteins % 1000:
                    print(f"Comparing {num_proteins} {pr.accession_number}")

                num_proteins += 1

                self._compare(pr.accession_number, project_original, project_process)

    def _compare(self, accession_number, project_original, project_process):
        protein_original = Protein.objects.get(
            project=project_original,
            accession_number = accession_number
        )

        try:
            protein_process = Protein.objects.get(
                project=project_process,
                accession_number = accession_number
            )
        except Exception as e:
            # TODO - why would this happen?
            print(f"Can't get process protein {accession_number}")
            return

        # Don't compare protein medians as they're not stored in the ICR json
        self._compare_protein_stat(PROTEIN_ABUNDANCES_RAW, protein_original, protein_process)
        self._compare_protein_stat(PROTEIN_ABUNDANCES_NORMALISED_MEDIAN, protein_original, protein_process)
        self._compare_protein_stat(PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN, protein_original, protein_process)
        self._compare_protein_stat(PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX, protein_original, protein_process)
        self._compare_protein_stat(PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX, protein_original, protein_process)
        self._compare_protein_stat(PROTEIN_ABUNDANCES_IMPUTED, protein_original, protein_process)
        self._compare_metrics(PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN, protein_original, protein_process, True)
        self._compare_metrics(PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX, protein_original, protein_process)



    def _compare_metrics(
        self,
        statistic_type_name,
        protein_original,
        protein_process,
        with_anova = False,
        dps = 0
    ):
        _, stat_prot_original = self._fetch_stats_type_and_stats(
            statistic_type_name,
            protein=protein_original
        )
        _, stat_prot_process = self._fetch_stats_type_and_stats(
            statistic_type_name,
            protein=protein_process
        )

        metrics_original = stat_prot_original.metrics
        metrics_process = stat_prot_process.metrics

        if not metrics_process:
            print(f"NO PROCESS METRICS FOR {statistic_type_name} for {protein_original.accession_number}")
            return

        for field in METRICS_STRING_FIELDS:
            if not metrics_process.get(field):
                print(f"No reading for METRICS for {statistic_type_name} for {protein_original.accession_number} for {field}")
                continue

            if metrics_original[field] != metrics_process[field]:
                print(f"NO MATCH {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")
            # else:
            #     print(f"Match for {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")

        for field in METRICS_NUMBER_FIELDS:
            if not metrics_process.get(field):
                print(f"No reading for METRICS for {statistic_type_name} for {protein_original.accession_number} for {field}")
                continue

            if round(metrics_original[field], dps) != round(metrics_process[field], dps):
                print(f"NO METRICS MATCH {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")
            # else:
            #     print(f"Metrics match for {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")

        if with_anova:
            if not metrics_process.get(ANOVA):
                print(f"NO PROCESS METRICS FOR {statistic_type_name} for {protein_original.accession_number}")
            else:
                for field in METRICS_ANOVA_FIELDS:
                    if not metrics_process[ANOVA].get(field):
                        print(f"No ANOVA reading for METRICS for {statistic_type_name} for {protein_original.accession_number} for {field}")
                        continue

                    if round(metrics_original[ANOVA][field], dps) != round(metrics_process[ANOVA][field], dps):
                        print(f"NO ANOVA METRICS MATCH {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[ANOVA][field]} vs {metrics_process[ANOVA][field]}")
                    # else:
                    #     print(f"ANOVA metrics match for {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[ANOVA][field]} vs {metrics_process[ANOVA][field]}")



    def _compare_protein_stat(
        self,
        statistic_type_name,
        protein_original,
        protein_process,
        dps = 0
    ):
        _, stat_prot_original = self._fetch_stats_type_and_stats(
            statistic_type_name,
            protein=protein_original
        )
        _, stat_prot_process = self._fetch_stats_type_and_stats(
            statistic_type_name,
            protein=protein_process
        )

        stat_prot_abundances_original = Abundance.objects.filter(
            statistic = stat_prot_original
        )
        stat_prot_abundances_process = Abundance.objects.filter(
            statistic = stat_prot_process
        )

        if len(stat_prot_abundances_original) != len(stat_prot_abundances_process):
            print(f"Protein abundance lengths don't match for {protein_original.accession_number}: {len(stat_prot_abundances_original)} vs {len(stat_prot_abundances_process)}")

        for ab_original in stat_prot_abundances_original:
            ab_process = stat_prot_abundances_process.filter(
                statistic = stat_prot_process,
                replicate__name = ab_original.replicate.name,
                sample_stage__name = ab_original.sample_stage.name
            ).first()

            if not ab_process:
                print(f"No reading for {statistic_type_name} for {protein_original.accession_number} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {ab_original.reading}")
                continue

            if round(ab_original.reading, dps) != round(ab_process.reading, dps):
                print(f"NO MATCH {statistic_type_name} for {protein_original.accession_number} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {round(ab_original.reading, dps)} vs {round(ab_process.reading, dps)}")
            # else:
            #     print(f"Match for {statistic_type_name} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {round(ab_original.reading, 1)} vs {round(ab_process.reading, 1)}")


    # TODO - copied from import_original
    def _fetch_stats_type_and_stats(self, statistic_type_name, project = None, protein = None, phospho = None):
        statistic_type = StatisticType.objects.get(name=statistic_type_name)

        if project:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, project=project)
        elif protein:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, protein=protein)
        elif phospho:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, phospho=phospho)
        else:
            raise Exception(f"_clear_and_fetch_stats needs a project, protein or phospho")

        return statistic_type, stat
