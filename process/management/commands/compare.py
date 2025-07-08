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
        self._compare_numbers(ABUNDANCES_RAW, protein_original, protein_process, None, None)
        self._compare_numbers(ABUNDANCES_NORMALISED_MEDIAN, protein_original, protein_process, None, None)
        self._compare_numbers(ABUNDANCES_NORMALISED_LOG2_MEAN, protein_original, protein_process, None, None)
        self._compare_numbers(ABUNDANCES_NORMALISED_MIN_MAX, protein_original, protein_process, None, None)
        self._compare_numbers(ABUNDANCES_NORMALISED_ZERO_MAX, protein_original, protein_process, None, None)
        self._compare_numbers(ABUNDANCES_IMPUTED, protein_original, protein_process, None, None)
        self._compare_metrics(ABUNDANCES_NORMALISED_LOG2_MEAN, protein_original, protein_process, None, None, True)
        self._compare_metrics(ABUNDANCES_NORMALISED_ZERO_MAX, protein_original, protein_process, None, None)

        # Compare phosphos
        phosphos_original = Phospho.objects.filter(
            protein = protein_original
        )

        for phospho_original in phosphos_original:
            try:
                phospho_process = Phospho.objects.get(
                    protein__project = project_process,
                    protein__accession_number = accession_number,
                    mod = phospho_original.mod
                )
            except Exception as e:
                # TODO - why would this happen?
                print(f"Can't get process phospho {accession_number} {phospho_original.mod}")
                print(e)
                return

            # Don't compare phospho medians as they're not stored in the ICR json
            self._compare_numbers(ABUNDANCES_RAW, None, None, phospho_original, phospho_process)
            self._compare_numbers(ABUNDANCES_NORMALISED_MEDIAN, None, None, phospho_original, phospho_process)
            self._compare_numbers(ABUNDANCES_NORMALISED_LOG2_MEAN, None, None, phospho_original, phospho_process)
            self._compare_numbers(ABUNDANCES_NORMALISED_MIN_MAX, None, None, phospho_original, phospho_process)
            self._compare_numbers(ABUNDANCES_NORMALISED_ZERO_MAX, None, None, phospho_original, phospho_process)
            self._compare_numbers(ABUNDANCES_IMPUTED, None, None, phospho_original, phospho_process)
            self._compare_metrics(ABUNDANCES_NORMALISED_LOG2_MEAN, None, None, phospho_original, phospho_process, True)
            self._compare_metrics(ABUNDANCES_NORMALISED_ZERO_MAX, None, None, phospho_original, phospho_process)


    def _compare_metrics(
        self,
        statistic_type_name,
        protein_original,
        protein_process,
        phospho_original,
        phospho_process,
        with_anova = False,
        dps = 2
    ):
        if protein_original and protein_process:
            _, stat_original = self._fetch_stats_type_and_stats(
                statistic_type_name,
                protein=protein_original
            )
            _, stat_process = self._fetch_stats_type_and_stats(
                statistic_type_name,
                protein=protein_process
            )
        elif phospho_original and phospho_process:
            _, stat_original = self._fetch_stats_type_and_stats(
                statistic_type_name,
                phospho = phospho_original
            )
            _, stat_process = self._fetch_stats_type_and_stats(
                statistic_type_name,
                phospho = phospho_process
            )
        else:
            raise Exception("compare_metrics requires either a pair of proteins or a pair of phosphos")

        metrics_original = stat_original.metrics
        metrics_process = stat_process.metrics

        if not metrics_process:
            if protein_original:
                print(f"NO PROCESS METRICS FOR {statistic_type_name} for {protein_original.accession_number}")
            else:
                print(f"NO PROCESS METRICS FOR {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod}")

            return

        for field in METRICS_STRING_FIELDS:
            if not metrics_process.get(field):
                if protein_original:
                    print(f"No reading string for METRICS for {statistic_type_name} for {protein_original.accession_number} for {field}")
                else:
                    print(f"No reading string for METRICS for {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod} for {field}")

                continue

            if metrics_original[field] != metrics_process[field]:
                if protein_original:
                    print(f"NO MATCH {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")
                else:
                    print(f"NO MATCH {statistic_type_name} for {phospho_original.protein.accession_number} mode {phospho_original.mod} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")

            # else:
            #     print(f"Match for {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")

        for field in METRICS_NUMBER_FIELDS:
            if metrics_process.get(field) is None:
                if protein_original:
                    print(f"No reading number for METRICS for {statistic_type_name} for {protein_original.accession_number} for {field}")
                else:
                    print(f"No reading number for METRICS for {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod} for {field}")

                continue

            if self._not_same(metrics_original[field], metrics_process[field], dps):
                if protein_original:
                    print(f"NO METRICS MATCH {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")
                else:
                    print(f"NO METRICS MATCH {statistic_type_name} for {phospho_original.protein.accession_number} mode {phospho_original.mod} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")
            # else:
            #     print(f"Metrics match for {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[field]} vs {metrics_process[field]}")

        if with_anova:
            if metrics_process.get(ANOVA) is None:
                if protein_original:
                    print(f"No ANOVA for {statistic_type_name} for {protein_original.accession_number} for {field}")
                else:
                    print(f"No ANOVA for {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod} for {field}")
            else:
                for field in METRICS_ANOVA_FIELDS:
                    if metrics_process[ANOVA].get(field) is None:
                        if protein_original:
                            print(f"No ANOVA METRICS for {statistic_type_name} for {protein_original.accession_number} for {field}")
                        else:
                            print(f"No ANOVA METRICS for {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod} for {field}")
                        continue

                    if self._not_same(metrics_original[ANOVA][field], metrics_process[ANOVA][field], dps):
                        if protein_original:
                            print(f"NO ANOVA METRICS MATCH {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[ANOVA][field]} vs {metrics_process[ANOVA][field]}")
                        else:
                            print(f"NO ANOVA METRICS MATCH {statistic_type_name} for {phospho_original.protein.accession_number} mode {phospho_original.mod} for {field} reading {metrics_original[ANOVA][field]} vs {metrics_process[ANOVA][field]}")
                    # else:
                    #     print(f"ANOVA metrics match for {statistic_type_name} for {protein_original.accession_number} for {field} reading {metrics_original[ANOVA][field]} vs {metrics_process[ANOVA][field]}")


    def _not_same(self, num1, num2, dps):
        return (round(num1, dps) != round(num2, dps)) and (abs(num1 - num2) > TOLERANCE * abs(num1))

    def _compare_numbers(
        self,
        statistic_type_name,
        protein_original,
        protein_process,
        phospho_original,
        phospho_process,
        dps = 0
    ):
        
        if protein_original and protein_process:
            _, stat_original = self._fetch_stats_type_and_stats(
                statistic_type_name,
                protein=protein_original
            )
            _, stat_process = self._fetch_stats_type_and_stats(
                statistic_type_name,
                protein=protein_process
            )
        elif phospho_original and phospho_process:
            _, stat_original = self._fetch_stats_type_and_stats(
                statistic_type_name,
                phospho=phospho_original
            )
            _, stat_process = self._fetch_stats_type_and_stats(
                statistic_type_name,
                phospho=phospho_process
            )
        else:
            raise Exception("compare_numbers requires either a pair of proteins or a pair of phosphos")

        stat_abundances_original = Abundance.objects.filter(
            statistic = stat_original
        )
        stat_abundances_process = Abundance.objects.filter(
            statistic = stat_process
        )

        if len(stat_abundances_original) != len(stat_abundances_process):
            if protein_original:
                print(f"Protein abundance lengths don't match for {protein_original.accession_number}: {len(stat_abundances_original)} vs {len(stat_abundances_process)}")
            else:
                print(f"Protein abundance lengths don't match for {phospho_original.protein.accession_number}: {len(stat_abundances_original)} vs {len(stat_abundances_process)}")

        for ab_original in stat_abundances_original:
            ab_process = stat_abundances_process.filter(
                statistic = stat_process,
                replicate__name = ab_original.replicate.name,
                sample_stage__name = ab_original.sample_stage.name
            ).first()

            if not ab_process:
                if protein_original:
                    print(f"No reading for {statistic_type_name} for {protein_original.accession_number} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {ab_original.reading}")
                else:
                    print(f"No reading for {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {ab_original.reading}")

                continue

            if self._not_same(ab_original.reading, ab_process.reading, dps):
                if protein_original:
                    print(f"NO MATCH {statistic_type_name} for {protein_original.accession_number} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {ab_original.reading} vs {ab_process.reading}")
                else:
                    print(f"NO MATCH {statistic_type_name} for {phospho_original.protein.accession_number} mod {phospho_original.mod} for {ab_original.replicate.name} for {ab_original.sample_stage.name} reading {ab_original.reading} vs {ab_process.reading}")
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
