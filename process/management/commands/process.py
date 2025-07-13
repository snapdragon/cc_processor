import re
import json
import logging
import math
import statistics
import copy
import requests
from scipy import stats

import numpy as np
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from django.db.models.query import QuerySet
from scipy.stats import f_oneway, moment
from sklearn.metrics import r2_score
from sklearn import linear_model
from scipy import stats

# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector

from process.models import (
    ColumnName,
    Project,
    Phospho,
    Protein,
    Replicate,
    SampleStage,
    StatisticType,
    Statistic,
    Abundance
)

from process.constants import (
    RAW,
    METRICS,
    LOG2_MEAN,
    ZERO_MAX,
    ANOVA,
    Q_VALUE,
    FISHER_G,
    PROTEIN_ABUNDANCES,
    ABUNDANCE_AVERAGE,
    NORMALISED,
    MEDIAN,
    MIN_MAX,
    LOG2_ARREST,
    IMPUTED,
    P_VALUE,
    F_STATISTICS,
    PROTEIN_OSCILLATION_ABUNDANCES,
    POSITION_ABUNDANCES,
    PHOSPHORYLATION_ABUNDANCES,
    PHOSPHO_REGRESSION,
    PHOSPHORYLATION_SITE,
    CURVE_FOLD_CHANGE,
    PROTEIN_PHOSPHO_CORRELATION,
    PHOSPHO_PROTEIN_CFC_RATIO,
    G_STATISTIC,
    FREQUENCY,
    TOTAL_PROTEIN_INDEX_FILE,
    KINASE_PREDICTION,
    PEPTIDE_SEQ,
    PHOSPHO,
    PEPTIDE,
    PEPTOOLS_ANNOTATIONS,
    CONSENSUS_MOTIF_MATCH,
    KINASE_MOTIF_MATCH,
    PROTEIN_INFO,
    DEFAULT_FISHER_STATS,
    PROTEIN_INFO_FIELDS,
    GENE_NAME,
    PROTEIN_NAME,
    ABUNDANCES_RAW,
    PROTEIN_MEDIAN,
    PHOSPHO_MEDIAN,

    ABUNDANCES_RAW,
    ABUNDANCES_IMPUTED,

    ABUNDANCES_NORMALISED_ZERO_MAX,
    ABUNDANCES_NORMALISED_MEDIAN,
    ABUNDANCES_NORMALISED_MIN_MAX,
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ABUNDANCES_NORMALISED_LOG2_ARREST,

    PHOSPHO_REGRESSION_ZERO_MAX,
    PHOSPHO_REGRESSION_LOG2_MEAN,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_RAW,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_IMPUTED,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_ZERO_MAX,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_MEDIAN,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_MIN_MAX,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_ARREST,

    PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX,
    PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,

    ICR_ABUNDANCE_REP_1,
    ICR_ABUNDANCE_REP_2,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

# TOOD - what is this file? Where did it come from?
# TODO - whatever it is, a lot of it could be trimmed, and maybe put in the DB
with open(f"./data/{TOTAL_PROTEIN_INDEX_FILE}") as outfile:
    index_protein_names = json.load(outfile)

class Command(BaseCommand):
    help = "Processes all proteins and phosphoproteins for a given project"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to process",
        )
        parser.add_argument(
            "--accession-number",
            help="Optional accession number",
        )
        parser.add_argument(
            "--with-bugs",
            help="Run with the original ICR bugs",
            action="store_true",
        )
        parser.add_argument(
            "--calculate-protein-medians",
            help="Calculate and store protein medians",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-proteins",
            help="Calculate and store protein results",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-phosphos",
            help="Calculate and store phospho results",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-phospho-medians",
            help="Calculate and store phospho medians",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-batch",
            help="Calculate and store batch values",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-all",
            help="Calculate and store all values",
            action="store_true"
        )

    def handle(self, *args, **options):
        project_name = options["project"]
        accession_number = options["accession_number"]
        with_bugs = options["with_bugs"]
        calculate_protein_medians = options["calculate_protein_medians"]
        calculate_proteins = options["calculate_proteins"]
        calculate_phospho_medians = options["calculate_phospho_medians"]
        calculate_phosphos = options["calculate_phosphos"]
        calculate_batch = options["calculate_batch"]
        calculate_all = options["calculate_all"]

        if not project_name:
            raise Exception("You must provide a project name.")

        if with_bugs and project_name != "ICR":
            raise CommandError("Only an ICR project can run --with-bugs")

        logger.info(f"Processing for project {project_name}, with bugs {with_bugs}")

        # TODO - return if no calculate flags are set?

        project = Project.objects.get(name=project_name)
        replicates = Replicate.objects.filter(project=project, mean=False)
        sample_stages = SampleStage.objects.filter(project=project).order_by('rank')

        project.with_bugs = with_bugs
        project.save()

        if calculate_protein_medians or calculate_all:
            self._calculate_medians(
                project,
                replicates,
                sample_stages,
                True,
                with_bugs
            )

        # N.B. there must be protein medians for this to run
        if calculate_proteins or calculate_all:
            if accession_number:
                logger.info(f"Processing protein {accession_number}")

                protein = Protein.objects.get(
                    project=project,
                    accession_number=accession_number
                )

                self._calculate_abundances_metrics(
                    replicates,
                    sample_stages,
                    protein = protein,
                    phospho = None,
                    with_bugs = with_bugs,
                )
            else:
                logger.info("Processing all proteins")

                proteins = Protein.objects.filter(project=project).iterator(chunk_size=100)

                for i, protein in enumerate(proteins):
                    if not i % 1000:
                        logger.info(f"Calculating protein {i} {protein.accession_number}")

                    self._calculate_abundances_metrics(
                        replicates,
                        sample_stages,
                        protein = protein,
                        phospho = None,
                        with_bugs = with_bugs,
                    )

        if calculate_phospho_medians or calculate_all:
            self._calculate_medians(
                project,
                replicates,
                sample_stages,
                False,
                with_bugs
            )

        if calculate_phosphos or calculate_all:
            if accession_number:
                logger.info(f"Processing phospho protein {accession_number}")

                protein = Protein.objects.get(
                    project=project,
                    accession_number=accession_number
                )

                self._calculate_phosphos(
                    replicates,
                    sample_stages,
                    protein,
                    with_bugs,
                )

            else:
                logger.info("Processing all phosphos")

                proteins = Protein.objects.filter(project=project).iterator(chunk_size=100)

                for i, protein in enumerate(proteins):
                    if not i % 1000:
                        logger.info(f"Calculating phospho protein {i} {protein.accession_number}")

                    self._calculate_phosphos(
                        replicates,
                        sample_stages,
                        protein,
                        with_bugs,
                    )

        if calculate_batch or calculate_all:
        #     self._add_protein_annotations(run)

            self._calculate_protein_q_and_fisher_g(project, replicates, sample_stages)

            self._calculate_phospho_q_and_fisher_g(project, replicates, sample_stages)

            self._add_oscillations(project, replicates, sample_stages, with_bugs)

        #     self._add_phospho_regression(project, replicates, sample_stages, with_bugs)

            # TODO - figure out how, or if at all, to get this working for SL
            #   Actually it no longer works for ICR either, related to
            #   P04264 somehow
            # self._generate_kinase_predictions(run)




    def _calculate_phosphos(self, replicates, sample_stages, protein, with_bugs):
        phosphos = Phospho.objects.filter(protein = protein)

        for phospho in phosphos:
            self._calculate_abundances_metrics(
                replicates,
                sample_stages,
                None,
                phospho,
                with_bugs
            )

    # TODO - similar to other functionality, consolidate
    def _calculate_phospho_q_and_fisher_g(self, project, replicates, sample_stages):
        logger.info("Calculating phospho q and fisher G values")
        # TODO - blank all q_values in DB?

        fisher_g_stats = self.calculate_fisher_g(project, replicates, sample_stages, phospho = True)

        # Get all phospho abundances for this project for log2 mean
        statistics = Statistic.objects.filter(
            statistic_type__name = ABUNDANCES_NORMALISED_LOG2_MEAN,
            phospho__protein__project = project
        )

        anova_stats: dict = {}

        for statistic in statistics:
            site_key = self._get_site_key(statistic)

            anova_stats[site_key] = {
                P_VALUE: statistic.metrics[ANOVA][P_VALUE]
            }

        anova_stats = self._add_q_value(anova_stats)

        for statistic in statistics:
            # Determine Fisher G
            site_key = self._get_site_key(statistic)

            fisher_output = DEFAULT_FISHER_STATS

            if site_key in fisher_g_stats:
                fisher_output = fisher_g_stats[site_key]

            statistic.metrics[FISHER_G] = fisher_output

            # Determine ANOVA q value
            q_value = 1

            if anova_stats.get(site_key):
                q_value = anova_stats[site_key][Q_VALUE]

            statistic.metrics[ANOVA][Q_VALUE] = q_value

            statistic.save()



    def _calculate_protein_q_and_fisher_g(self, project, replicates, sample_stages):
        logger.info("Calculating q and fisher G values for proteins.")

        # TODO - blank all q_values in DB?

        # Get all protein abundances for this project for log2 mean
        statistics = Statistic.objects.filter(
            statistic_type__name = ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein__project = project
        )

        fisher_g_stats = self.calculate_fisher_g(project, replicates, sample_stages)

        anova_stats: dict = {}

        for statistic in statistics:
            anova_stats[statistic.protein] = {}
            anova_stats[statistic.protein][P_VALUE] = statistic.metrics[ANOVA][P_VALUE]

        anova_stats = self._add_q_value(anova_stats)

        for statistic in statistics:
            # Determine q value
            q_value = 1

            if anova_stats.get(statistic.protein):
                q_value = anova_stats[statistic.protein][Q_VALUE]

            statistic.metrics[ANOVA][Q_VALUE] = q_value

            # Determine Fisher G value
            fisher_output = DEFAULT_FISHER_STATS

            if statistic.protein.accession_number in fisher_g_stats:
                fisher_output = fisher_g_stats[statistic.protein.accession_number]

            statistic.metrics[FISHER_G] = fisher_output

            statistic.save()



    # createAbundanceDf
    def _create_abundance_dataframe(self, project, replicates, sample_stages, phospho, phospho_ab, phospho_reg):
        abundance_table = {}

        # TODO - could combine these two loops using a function to generate the keys
        if phospho:
            # Determine which phospho readings to use
            statistic_type_name = ABUNDANCES_NORMALISED_LOG2_MEAN

            if phospho_ab:
                statistic_type_name = PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN
            elif phospho_reg:
                statistic_type_name = PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN

            # Get all log2 mean abundances for all chosen phospho readings for this project
            abundances = Abundance.objects.filter(
                statistic__statistic_type__name = statistic_type_name,
                statistic__phospho__protein__project = project,
            ).iterator(chunk_size=100)

            # abundances = Abundance.objects.filter(
            #     statistic__statistic_type__name = statistic_type_name,
            #     statistic__phospho__protein__project = project,
            # )[:10000]


            num_abundances = 0

            for abundance in abundances:
                if not num_abundances % 100000:
                    logger.info(f"Processing dataframe phospho abundance {num_abundances}")

                num_abundances += 1

                site_key = self._get_site_key(abundance.statistic)

                if abundance_table.get(site_key) is None:
                    abundance_table[site_key] = {}

                rep_timepoint = f"{abundance.replicate.name}_{abundance.sample_stage.name}"

                abundance_table[site_key][rep_timepoint] = abundance.reading
        else:
            # Get all log2 mean abundances for all proteins for this project
            abundances = Abundance.objects.filter(
                statistic__statistic_type__name = ABUNDANCES_NORMALISED_LOG2_MEAN,
                statistic__protein__project = project,
            ).iterator(chunk_size=100)

            num_abundances = 0

            for abundance in abundances:
                if not num_abundances % 100000:
                    logger.info(f"Processing dataframe protein abundance {num_abundances}")

                num_abundances += 1

                pan = abundance.statistic.protein.accession_number

                # Not all proteins have protein results, some are phospho only
                # TODO - is this necessary? Would phospho-only results have a record?
                if abundance.reading is None:
                    continue

                if abundance_table.get(pan) is None:
                    abundance_table[pan] = {}

                rep_stage_name = f"{abundance.replicate.name}_{abundance.sample_stage.name}"

                abundance_table[pan][rep_stage_name] = abundance.reading

        time_course_abundance_df = pd.DataFrame(abundance_table)
        time_course_abundance_df = time_course_abundance_df.T

        new_cols = []

        # TODO - is this ordering correct?
        for sample_stage in sample_stages:
            for replicate in replicates:
                new_cols.append(f"{replicate.name}_{sample_stage.name}")

        # Rearrange the column order so replicates are near eatch other
        try:
            time_course_abundance_df = time_course_abundance_df[new_cols]
        except Exception as e:
            logger.error("DF FAILED")
            logger.error(phospho)
            logger.error(new_cols)
            logger.error(time_course_abundance_df.columns)
            logger.error(time_course_abundance_df)
            logger.error(e)
            exit()

        return time_course_abundance_df

    # TODO - is this really needed any more?
    def _tp(self, abundances):
        """
        Creates a list of each abundance of a stage name across replicates
        """
        res = []

        for abundance in abundances:
            if abundance.reading is not None:
                res.append(abundance.reading)

        # for replicate in replicates:
        #     if stage_name in readings[replicate.name]:
        #         reading = readings[replicate.name][stage_name]

        #         if reading is not None:
        #             res.append(float(reading))

        return res
    
    # TODO - straight up lifted from ICR, simplify ideally using a library
    def _calculate_ANOVA(self, statistic_type_name, sample_stages, protein = None, phospho = None):
        # Defaults if not enough replicates
        p_value = 1
        f_statistic = 1

        stat_log2_mean = self._get_statistic(
            statistic_type_name,
            protein = protein,
            phospho = phospho,
        )

        abundances = Abundance.objects.filter(
            statistic = stat_log2_mean,
            replicate__mean = False
        ).order_by(
            'sample_stage__rank'
        )

        # TODO - not a good name
        stage_names = []

        try:
            for sample_stage in sample_stages:
                stage_names.append(
                    self._tp(
                        abundances.filter(sample_stage=sample_stage),
                    )
                )

            # Each entry must have at least two points for f_oneway to work    
            timepoints = [x for x in stage_names if x != [] and len(x) > 1]

            if len(timepoints) > 1:
                # TODO - study f_oneway
                one_way_anova = f_oneway(*timepoints)

                f_statistic = one_way_anova[0].item()
                p_value = one_way_anova[1].item()

                if np.isnan(p_value):
                    p_value = 1

            metrics = stat_log2_mean.metrics
            metrics[ANOVA] = {
                P_VALUE: p_value,
                F_STATISTICS: f_statistic
            }

            stat_log2_mean.save()

        except Exception as e:
            logger.error("Error in _calculate_ANOVA")
            logger.error(e)

    # TODO - lifted from ICR, change names
    def _calculate_metrics(
        self,
        statistic_type_name: str,
        replicates,
        sample_stages,
        protein = None,
        phospho = None,
    ):
        metrics = {}

        if protein:
            stat = self._get_statistic(statistic_type_name, protein = protein)
        else:
            stat = self._get_statistic(statistic_type_name, phospho = phospho)

        abundances_for_readings = Abundance.objects.filter(
            statistic=stat
        ).order_by(
            "replicate__id",
            "sample_stage__rank"
        )

        # Clear previous metrics
        stat.metrics = {}
        stat.save()

        # TODO - converts abundances to the old format used by this function.
        #   Refactor this function to use abundances directly.
        readings = {}
        readings_averages = {}

        for abfr in abundances_for_readings:
            if abfr.replicate.mean:
                readings_averages[abfr.sample_stage.name] = abfr.reading

            if readings.get(abfr.replicate.name) is None:
                readings[abfr.replicate.name] = {}

            readings[abfr.replicate.name][abfr.sample_stage.name] = abfr.reading


        # TODO - all old code below. Tidy this.        
        metrics = {}

        abundances = []
        # TODO - how does ICR cope with Nones but not this? Does it use imputed values?
        abundance_averages = {
            key:val for key, val in readings_averages.items() if val is not None
        }
        abundance_averages_list = [
            val for val in abundance_averages.values()
        ]

        for replicate in replicates:
            for sample_stage in sample_stages:
                # Not all replicates are populated. The qc process removes any replicates
                #   with sample stage values that are all None.
                if readings.get(replicate.name) and not replicate.mean:
                    if abundance := readings[replicate.name].get(sample_stage.name):
                        # TODO - sometimes the sample_stage.name isn't set. Why?
                        #   Is there something wrong with the populate script?
                        # TODO - why does this add all reps for standard deviation calculation?
                        #   Why isn't it on a per-replicate basis?
                        abundances.append(abundance)

        std = None

        if len(abundances) > 1:
            std = statistics.stdev(abundances)

        # TODO - is this the right thing to do? Go through all of this function
        #   to check the behaviour of its various parts.
        #   Almost all of the warnings output are from this function.
        if not len(abundance_averages_list):
            return

        try:
            residuals_array = np.polyfit(
                range(0, len(abundance_averages)),
                abundance_averages_list,
                2,
                full=True,
            )[1]

            if len(residuals_array) == 0:
                # TODO - what? And get rid of comment below.
                # eg Q9HBL0 {'G2_2': 0.4496, 'G2/M_1': 0.7425, 'M/Early G1': 1.0}
                residuals = 5
            else:
                residuals = np.polyfit(
                    range(0, len(abundance_averages)),
                    abundance_averages_list,
                    2,
                    full=True,
                )[1][0]

            r_squared = self._polyfit(
                range(0, len(abundance_averages)), abundance_averages_list, 2
            )

            # TODO - why does this happen?
            # Converted to None as NaN can't be turned to json.
            if math.isnan(r_squared) or math.isinf(r_squared):
                r_squared = None

            max_fold_change = max(abundance_averages.values()) - min(
                abundance_averages.values()
            )

            metrics = {
                "variance": moment(abundance_averages_list, moment=2),
                "skewness": moment(abundance_averages_list, moment=3),
                "kurtosis": moment(abundance_averages_list, moment=4),
                "peak": max(abundance_averages, key=abundance_averages.get),
                "max_fold_change": max_fold_change,
                "residuals": residuals,
                "R_squared": r_squared,
            }

            # if we have info for at least 2 replicates
            # TODO - why does it get the values for the second one?
            if len(readings) >= 2:
                curve_fold_change, curve_peak = self._calculate_curve_fold_change(
                    readings, replicates
                )

                residuals_all, r_squared_all = self._calculate_residuals_R2_all(readings, replicates)

                metrics = {
                    "standard_deviation": std, 
                    "variance_average": round(moment(abundance_averages_list, moment=2),2),
                    "skewness_average": moment(abundance_averages_list, moment=3),
                    "kurtosis_average": moment(abundance_averages_list, moment=4),
                    "peak_average": max(abundance_averages, key=abundance_averages.get),
                    "max_fold_change_average": max_fold_change,
                    "residuals_average": residuals,
                    "R_squared_average": r_squared,
                    "residuals_all": residuals_all,
                    "R_squared_all": r_squared_all,
                    CURVE_FOLD_CHANGE: curve_fold_change,
                    "curve_peak": curve_peak,
                }

        except Exception as e:
            logger.error("Error in _calculate_metrics")
            logger.error(e)

        stat.metrics = metrics
        stat.save()


    # TODO - this is straight up lifted from ICR, change and ideally use a library
    # _calcResidualsR2All
    # TODO - tested
    def _calculate_residuals_R2_all(self, readings, replicates):
        residuals_all = None
        r_squared_all = None

        x, y, _ = self._generate_xs_ys(readings, replicates)

        if len(x) != len(y):
            return residuals_all, r_squared_all

        p = np.poly1d(np.polyfit(x, y, 2))
        curve_abundances = p(x)
        residuals_all = np.polyfit(x, y, 2, full=True)[1][0]
        r_squared_all = round(r2_score(y, curve_abundances), 2)

        return residuals_all.item(), r_squared_all

    # TODO - this is straight up lifted from ICR, change and ideally use a library
    # TODO - tested (ish)
    def _calculate_curve_fold_change(self, readings, replicates):
        """
        Calculates the curve_fold_change and curve peaks for the three or two replicates normalised abundance for each protein.
        """
        curve_fold_change = None
        curve_peak = None

        x, y, stage_names_map = self._generate_xs_ys(readings, replicates)

        if len(x) == len(y):
            p = np.poly1d(np.polyfit(x, y, 2))
            curve_abundances = p(x)

            # find the timepoint peak of the curve
            curve_index = x[list(curve_abundances).index(max(curve_abundances))]
            for time_point, index in stage_names_map.items():
                if index == curve_index:
                    curve_peak = time_point

            # Calculate the fold change from the curve
            curve_fold_change = max(curve_abundances) / max(
                0.05, min(curve_abundances)
            )
            curve_fold_change = curve_fold_change.item()

        return curve_fold_change, curve_peak

    # TODO - this is straight up lifted from ICR. Replace it, ideally with a library call
    # TODO - write a test for it first
    def _polyfit(self, x, y, degree):
        coeffs = np.polyfit(x, y, degree)
        p = np.poly1d(coeffs)
        yhat = p(x)
        ybar = np.mean(y)
        ssres = np.sum((y - yhat) ** 2)
        sstot = np.sum((y - ybar) ** 2)
        r_squared = 1 - (ssres / sstot)
        return round(r_squared, 2)

    # TODO - what is this for? Why does it work on ABUNDANCES_NORMALISED_MIN_MAX?
    def _impute(
        # TODO - all these pr types are wrong, and also probably bad variable names
        self,
        replicates: QuerySet[Replicate],
        protein = None,
        phospho = None,
    ):
        _, stat_imputed = self._clear_and_fetch_stats(
            ABUNDANCES_IMPUTED,
            protein = protein,
            phospho = phospho
        )

        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MIN_MAX,
            protein = protein,
            phospho = phospho,
        ).order_by(
            "sample_stage__rank"
        )

        abundances_by_rep = {}

        for abundance in abundances:
            if abundances_by_rep.get(abundance.replicate) is None:
                abundances_by_rep[abundance.replicate] = []

            abundances_by_rep[abundance.replicate].append(abundance)

        for replicate in abundances_by_rep:
            abs = abundances_by_rep[replicate]

            for i, abundance in enumerate(abs):
                reading = 0

                if abundance.reading is not None:
                    reading = abundance.reading
                else:
                    last = None
                    next = None

                    # Go backwards to find most recent non-None value
                    for j in range(i - 1, -1, -1):
                        prev_abundance = abs[j]

                        if prev_abundance.reading is not None:
                            last = (i - j, prev_abundance.reading)
                            break

                    # Go forward to find next value
                    for j in range(i + 1, len(abs)):
                        # TODO - is it right that it loops back to the beginning?
                        #   Why doesn't going backwards loop too?
                        next_abundance = abs[j % len(abs)]

                        if next_abundance.reading is not None:
                            next = (j, next_abundance.reading)
                            break

                    if last and next:
                        # Linear imputation between nearest timepoints
                        last_offset, last_reading = last
                        next_offset, next_reading = next
                        step_height = (last_reading - next_reading) / (last_offset + next_offset)
                        reading = next_offset * step_height + next_reading

                Abundance.objects.create(
                    statistic=stat_imputed,
                    replicate=abundance.replicate,
                    sample_stage=abundance.sample_stage,
                    reading=reading
                )

    def _calculate_zero_or_min_normalisation(
        self,
        replicates,
        protein = None,
        phospho = None,
        zero_min = False
    ):
        # TODO - can be tidied, only constants are different
        if zero_min:
            _, stat = self._clear_and_fetch_stats(
                ABUNDANCES_NORMALISED_ZERO_MAX,
                protein = protein,
                phospho = phospho
            )

            abundances = self._get_abundances(
                ABUNDANCES_NORMALISED_MEDIAN,
                protein = protein,
                phospho = phospho
            )
        else:
            _, stat = self._clear_and_fetch_stats(
                ABUNDANCES_NORMALISED_MIN_MAX,
                protein = protein,
                phospho = phospho
            )

            abundances = self._get_abundances(
                ABUNDANCES_NORMALISED_LOG2_MEAN,
                protein = protein,
                phospho = phospho
            )
    
        for replicate in replicates:
            abs = abundances.filter(replicate=replicate)

            min_value = 0
            max_value = 0

            readings = []

            for ab in abs:
                readings.append(ab.reading)

            if len(readings):
                if not zero_min:
                    min_value = min(readings)

                max_value = max(readings)
            
            denominator = max_value - min_value

            for ab in abs:
                reading = ab.reading

                if reading is None or denominator == 0:
                    # TODO - why 0.5?
                    reading_normalised = 0.5
                else:
                    # TODO - why rounded?
                    reading_normalised = round((reading - min_value) / denominator, 4)

                Abundance.objects.create(
                    statistic=stat,
                    replicate=ab.replicate,
                    sample_stage=ab.sample_stage,
                    reading=reading_normalised
                )


    # calclog2RelativeAbundance
    def _calculate_relative_log2_normalisation(self, protein = None, phospho = None):
        _, stat_normalised_log2_mean = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein = protein,
            phospho = phospho
        )


        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MEDIAN,
            protein = protein,
            phospho = phospho,
        )

        total_abundances = 0
        total_lengths = 0

        for abundance in abundances:
            if abundance.reading is not None:
                total_abundances += math.log2(abundance.reading)
                total_lengths += 1

        mean = None

        if total_lengths != 0:
            mean = total_abundances / total_lengths

        for abundance in abundances:
            normalised_abundance = self._round(math.log2(abundance.reading) - mean)

            Abundance.objects.create(
                statistic=stat_normalised_log2_mean,
                replicate=abundance.replicate,
                sample_stage=abundance.sample_stage,
                reading=normalised_abundance
            )



    def _calculate_arrest_log2_normalisation(self, protein = None, phospho = None):
        ARRESTING_AGENT = "Nocodozole"

        _, stat_normalised_log2_arrest = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_LOG2_ARREST,
            protein = protein,
            phospho = phospho
        )

        if protein:
            # TODO - this is a hack, add the field to the Project model.
            if protein.project.name == "ICR":
                ARRESTING_AGENT = "Palbo"
        else:
            if phospho.protein.project.name == "ICR":
                ARRESTING_AGENT = "Palbo"

        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MEDIAN,
            protein = protein,
            phospho = phospho
        )

        for abundance in abundances:
            reading = abundance.reading

            # Get the arresting agent abundance for this replicate
            arrest_abundance = abundances.filter(
                sample_stage__name=ARRESTING_AGENT,
                replicate=abundance.replicate
            ).first()

            # Not all replicates have an arresting agent value
            # TODO - do we need to check for the other two Nones?
            if reading is not None and arrest_abundance is not None and arrest_abundance.reading is not None:
                log2_reading = self._round(math.log2(
                    reading / arrest_abundance.reading
                ))

                Abundance.objects.create(
                    statistic=stat_normalised_log2_arrest,
                    replicate=abundance.replicate,
                    sample_stage=abundance.sample_stage,
                    reading=log2_reading
                )


    def _get_statistic_type(self, name):
        return StatisticType.objects.get(name=name)

    def _get_statistic(self, statistic_type_name, project=None, protein=None, phospho=None):
        statistic_type = self._get_statistic_type(statistic_type_name)

        if project:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, project=project)
        elif protein:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, protein=protein)
        elif phospho:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, phospho=phospho)
        else:
            raise Exception(f"_get_protein_medians needs a project, protein or phospho")

        return stat

    def _get_protein_readings(self, protein):
        return self._get_abundances(ABUNDANCES_RAW, protein=protein)

    def _get_protein_medians(self, project):
        return self._get_abundances(PROTEIN_MEDIAN, project=project)

    def _get_abundances(self, statistic_type_name, project=None, protein=None, phospho=None):
        statistic = self._get_statistic(statistic_type_name, project, protein, phospho)

        return Abundance.objects.filter(statistic=statistic)

    def _get_abundance(self, statistic, replicate, sample_stage):
        return Abundance.objects.get(statistic=statistic, replicate=replicate, sample_stage=sample_stage)

    # TODO - only used once, remove
    def _delete_abundances(self, statistic):
        Abundance.objects.filter(statistic=statistic).delete()

    def _calculate_normalised_medians(self, protein = None, phospho = None):
        _, stat_normalised_medians = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_MEDIAN,
            protein = protein,
            phospho = phospho
        )

        if protein:
            stat_medians = self._get_statistic(PROTEIN_MEDIAN, project = protein.project)
        else:
            stat_medians = self._get_statistic(PHOSPHO_MEDIAN, project = phospho.protein.project)

        readings = self._get_abundances(ABUNDANCES_RAW, protein=protein, phospho=phospho)

        # TODO - get rid of _get_protein_readings? It may not be that useful

        for prr in readings:
            reading = prr.reading

            # TODO - is this if statement needed?
            if reading is not None:
                median = self._get_abundance(
                    stat_medians, 
                    prr.replicate,
                    prr.sample_stage
                )

                normalised_reading = reading / median.reading

                Abundance.objects.create(
                    statistic=stat_normalised_medians,
                    replicate=prr.replicate,
                    sample_stage=prr.sample_stage,
                    reading=normalised_reading
                )


    # calculateAverageRepAbundance
    def _calculate_means(
        self,
        statistic_type_name: str,
        protein = None,
        phospho = None,
        with_bugs: bool = False,
        imputed: bool = False,
    ):
        readings: dict = {}

        abundances = self._get_abundances(
            statistic_type_name,
            protein = protein,
            phospho = phospho
        )

        for abundance in abundances:
            sample_stage = abundance.sample_stage

            if readings.get(sample_stage) is None:
                readings[sample_stage] = []

            if with_bugs and not imputed:
                # We throw away the second reading
                # TODO - how will this behave towards None?
                if len(readings[sample_stage]) == 1:
                    continue

            if abundance.reading is not None:
                readings[sample_stage].append(abundance.reading)

        if protein:
            stat = self._get_statistic(statistic_type_name, protein=protein)
        else:
            stat = self._get_statistic(statistic_type_name, phospho=phospho)

        # Delete all mean abundances for this statistic
        Abundance.objects.filter(statistic=stat, replicate__mean=True).delete()

        if protein:
            mean_replicate = Replicate.objects.get(project=protein.project, mean=True)
        else:
            mean_replicate = Replicate.objects.get(project=phospho.protein.project, mean=True)

        for sample_stage in readings:
            reading_list = readings[sample_stage]

            if len(reading_list):
                mean = sum(reading_list) / len(reading_list)
                mean = self._round(mean)
            else:
                # TODO - is this the right thing to do?
                mean = None

            Abundance.objects.create(
                statistic=stat,
                replicate=mean_replicate,
                sample_stage=sample_stage,
                reading=mean
            )


    # TODO - think this could be tidied
    def _generate_xs_ys(self, readings, replicates):
        # TODO - why do we use the second replicate? It's in ICR
        second_replicate = list(replicates)[1]

        stage_names_map = {}

        # TODO - change this to use samples_stages list
        for i, stage_name in enumerate(readings[second_replicate.name].keys()):
            stage_names_map[stage_name] = i

        x = []
        for stage_name in readings.get(second_replicate.name, {}):
            x.append(stage_names_map[stage_name])
        x.sort()

        y = []
        for stage_name in stage_names_map:
            value = readings.get(second_replicate.name, {}).get(stage_name)
            if value is not None:
                y.append(value)

        return x, y, stage_names_map
    
    # TODO - lifted from ICR, change
    def _get_consensus_kinase_pred(self, uniprot_accession, phosphosite):
        """
        Predicts the kinase most likely to phosphorylate a phosphorylation site 
        based on the consensus approach.
        """
        phospho_kinases_class = {}
        peptide_seq = self._get_peptide_aligned(uniprot_accession, phosphosite)

        motifs = [
            {"motif": "plk", "pattern": ".{3}[DNE].[ST][FGAVLIMW][GAVLIPFMW].{3}"},
            {"motif": "cdk", "pattern": ".{3}..[ST]P.[RK]."},
            {"motif": "aurora", "pattern": ".{3}R.[ST][GAVLIFMW].{4}"},
            {"motif": "stp_consensus", "pattern": ".{5}[ST]P.{4}"},
            {"motif": "stq_consensus", "pattern": ".{5}[ST]Q.{4}"},
            {"motif": "krxst_consensus", "pattern": ".{3}[KR].[ST].{5}"},
            {"motif": "krxxst_consensus", "pattern": ".{2}[KR].{2}[ST].{5}"},
        ]

        total_not_matched = 0
        matches = 0
        motif_matches = []
        for m in motifs:
            motif = m["motif"]
            pattern = m["pattern"]
            res = re.match(pattern, peptide_seq)
            if res:
                matches+=1
                motif_matches.append(motif)
            # TODO - rename to accession_number?
            phospho_kinases_class = {"accession":uniprot_accession,"site":phosphosite, "peptide_seq": peptide_seq, "kinase_motif_match": motif_matches}
        if matches == 0:
            total_not_matched+=1
            phospho_kinases_class = {"accession":uniprot_accession, "site":phosphosite, "peptide_seq": peptide_seq, "kinase_motif_match": ["-"]}

        return phospho_kinases_class

    # TODO - lifted from ICR
    # getPeptideAligned
    def _get_peptide_aligned(self, uniprot_accession, phosphosite):
        """
        Checks if the phosphosite is centered in the peptide sequence and it aligns it if not.
        """
        aa = phosphosite[0]
        position = int(phosphosite[1::]) -1
        site =  str(position+1)

        # TODO - this was put in as without it the line
        #       if not len(peptide_seq):
        #   failed due to lack of access to the variable. This bug is in the
        #   original, how did it ever work? Was it never called?
        #   Why did it suddenly start erroring? It happened just after I
        #   first removed some of the R code from calculate_fisher_g.
        #   It's something to do with P04264, one of the duplicated proteins
        #   for ICR. It's the reason I made it like get_or_create in import_protein
        #   on duplicate, not just continue.
        # peptide_seq = []

        if PHOSPHO in index_protein_names[uniprot_accession]:
            if site in index_protein_names[uniprot_accession][PHOSPHO]:
                if PEPTIDE_SEQ in index_protein_names[uniprot_accession][PHOSPHO][site]:
                    peptide_seq = index_protein_names[uniprot_accession][PHOSPHO][site][PEPTIDE_SEQ]
                else:
                    peptide_seq = self._get_peptide_sequence(uniprot_accession, phosphosite)

        phospho_alignment = ""

        #Discard invalid inputs
        if not len(peptide_seq):
            peptide_seq = index_protein_names[uniprot_accession][PHOSPHO][site][PEPTIDE]
            phospho_alignment = peptide_seq[5:16]

        # Middle of the protein sequence
        elif len(peptide_seq) == 11:
            if aa == peptide_seq[5]:
                phospho_alignment = peptide_seq
            else:
                # Site not in the middle of seq
                peptide_seq_new = index_protein_names[uniprot_accession][PHOSPHO][site][PEPTIDE]
                phospho_alignment = peptide_seq_new[5:16]
                if  uniprot_accession == "P62861" and phosphosite == 'S5':
                    phospho_alignment = "-KVHGSLARAG"

        # Missing Positions
        elif len(peptide_seq) < 11:
            # beginning of the protein sequence
            if position < 5:
                increase = 5 - position
                peptide_seq = ("-" * increase) + peptide_seq
                phospho_alignment = peptide_seq
            else:
                # end of the protein sequence
                increase = 11 - len(peptide_seq)
                peptide_seq = peptide_seq + ("-" * increase)
                if aa != peptide_seq[5]:
                    peptide_seq = peptide_seq[:5] + aa + peptide_seq[5 + 1:]
                phospho_alignment = peptide_seq

        return phospho_alignment

    # TODO - lifted from ICR
    def _get_peptide_sequence(self, uniprot_accession, phosphosite):
        """
        Creates a peptide sequence which is a substring of the original protein sequence.
        +5/-5
        amino acids from the phosphorylation site.
        """
        sequence = self._get_protein_seq(uniprot_accession)

        # Start counting from 1
        position = int(phosphosite[1::]) -1
        sequence_len = len(sequence)
    
        start = 0
        end = sequence_len-1

        # start of protein sequence
        if position < 5:
            start = 0
            end = position + 6
        # middle of protein sequence
        elif sequence_len > position + 6:
            start = position - 5
            end = position + 6
        # end of protein sequence
        elif  sequence_len < position + 6 or sequence_len == position + 6:
            start = position - 5
            end = sequence_len

        peptide_sequence = sequence[start:end]

        return peptide_sequence

    # TODO - lifted from ICR
    def _get_protein_sequence(self, uniprot_accession):
        """
        Fetches the protein sequence for the given UniProt accession by querying UniProt directly.
        """
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.fasta"
    
        try:
            response = requests.get(url)
            if response.status_code == 200:
                lines = response.text.splitlines()
                sequence = ''.join(lines[1:]) 
                return sequence
            else:
                return f"Error: Unable to fetch data for accession {uniprot_accession}. HTTP Status: {response.status_code}"
        except requests.RequestException as e:
            return f"Error: Unable to connect to UniProt. {str(e)}"



    # addProteinOscillations
    def _add_oscillations(self, project, replicates, sample_stages, with_bugs):
        # TODO - check all logging statements for similarity to ICR
        logger.info("Adding Protein Oscillation Normalised Abundances")

        self._generate_protein_oscillation_metrics(project, replicates, sample_stages, with_bugs)

        fisher_g_stats = self.calculate_fisher_g(project, replicates, sample_stages, phospho = True, phospho_ab = True)

        ps_and_qs = {}

        # run_results = self._fetch_run_results(run)

        num_proteins = 0

        statistics = Statistic.objects.filter(
            statistic_type__name=PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
            phospho__protein__project=project
        ).iterator(chunk_size=100)

        # for rr in run_results:
        for statistic in statistics:
            if not num_proteins % 1000:
                print(f"Calculating protein oscillation {num_proteins} {statistic.phospho.protein.accession_number}")
 
            num_proteins += 1

            # for mod in pprpa:
            # TODO - generate with a site_key function?
            phospho_key = f"{statistic.phospho.protein.accession_number}_{statistic.phospho.phosphosite}"

            if phospho_key not in ps_and_qs:
                ps_and_qs[phospho_key] = {}

            # print("++++ STAT")
            # print(statistic)
            # print(statistic.metrics)

            if (statistic.metrics.get(ANOVA) is not None) and (statistic.metrics[ANOVA].get(P_VALUE) is not None):
                # print(f"Match {statistic.metrics[ANOVA][P_VALUE]}")

                ps_and_qs[phospho_key][P_VALUE] = statistic.metrics[ANOVA][P_VALUE]

            # pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            # # for mod in pprpa:
            # # TODO - generate with a site_key function?
            # phospho_key = f"{phospho.protein.accession_number}_{phospho.phosphosite}"

            # if phospho_key not in ps_and_qs:
            #     ps_and_qs[phospho_key] = {}

            # ps_and_qs[phospho_key][P_VALUE] = pprpa[mod][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][P_VALUE]

        ps_and_qs = self._add_q_value(ps_and_qs)

        # run_results = self._fetch_run_results(run)

        num_proteins = 0

        statistics = Statistic.objects.filter(
            statistic_type__name=PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
            phospho__protein__project=project
        ).iterator(chunk_size=100)

        # for rr in run_results:
        # TODO - this code is very similar to others, consolidate?
        for statistic in statistics:
            if not num_proteins % 1000:
                print(f"Calculating protein oscillation {num_proteins} {statistic.phospho.protein.accession_number}")
 
            num_proteins += 1

        # for rr in run_results:
        #     pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

        #     # TODO - why check for raw?
        #     if not len(pprpa) or not len(rr.combined_result[PROTEIN_ABUNDANCES][RAW]):
        #         continue

            # for mod in pprpa:
            #     if pprpa[mod].get(PROTEIN_OSCILLATION_ABUNDANCES) is None:
            #         continue

            # mod_key = f"{rr.protein.accession_number}_{pprpa[mod][PHOSPHORYLATION_SITE]}"
            phospho_key = f"{statistic.phospho.protein.accession_number}_{statistic.phospho.phosphosite}"

            q_value = 1

            if phospho_key in ps_and_qs:
                q_value = ps_and_qs[phospho_key][Q_VALUE]

            statistic.metrics[ANOVA][Q_VALUE] = q_value

            fisher_output = DEFAULT_FISHER_STATS

            if phospho_key in fisher_g_stats:
                fisher_output = fisher_g_stats[phospho_key]

            statistic.metrics[FISHER_G] = fisher_output

            statistic.save()

            #     pprpa[mod][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS]["Fisher_G"] = fisher_output

            # rr.save()



    # addPhosphoRegression
    def _add_phospho_regression(self, project, replicates, sample_stages, with_bugs):
        # TODO - check and revise all comments
        """
        Normalise the phospho abundance on the protein abundance
        Calculates and adds all the Regressed Phospho Abundance and their metrics for each phosphosite.
        Phospho = Dependent = Y
        Protein = Independent = X
        Linear Model => y = ax + b
        Residuals = Y - Y_predict
        """
        logger.info("Adding Phospho Normalised on Protein Abundances - Regression")

        self._generate_phospho_regression_metrics(project, replicates, sample_stages, with_bugs)

        # Fisher G Statistic - Phospho
        time_course_fisher_dict = self.calculate_fisher_g(project, replicates, sample_stages, phospho = True, phospho_ab = False, phospho_reg = True)

        regression_info = {}

        run_results = self._fetch_run_results(run)

        num_proteins = 0

        # TODO - looks very close to the code in protein oscillation. Make a function?
        for rr in run_results:
            if not num_proteins % 1000:
                print(f"Adding phospho regression {num_proteins} for {rr.protein.accession_number}")

            num_proteins += 1

            pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            if not len(pprpa):
                continue

            for mod in pprpa:
                if PHOSPHO_REGRESSION in pprpa[mod]:
                    # TODO - a hack
                    phospho_key = f"{rr.protein.accession_number}_{pprpa[mod][PHOSPHORYLATION_SITE]}"

                    if phospho_key not in regression_info:
                        regression_info[phospho_key] = {}

                    if pprpa[mod][PHOSPHO_REGRESSION][LOG2_MEAN].get(METRICS):
                        regression_info[phospho_key][P_VALUE] = pprpa[mod][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][P_VALUE]

        regression_info = self._add_q_value(regression_info)

        run_results = self._fetch_run_results(run)

        for rr in run_results:
            pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            if not len(pprpa) or not len(rr.combined_result[PROTEIN_ABUNDANCES][RAW]):
                continue

            # TODO - again, why called phosphosite and not mod?
            for phosphosite in pprpa:
                if PHOSPHO_REGRESSION in pprpa[phosphosite]:
                    mod_key = f"{rr.protein.accession_number}_{pprpa[phosphosite][PHOSPHORYLATION_SITE]}"

                    q_value = 1

                    if mod_key in regression_info:
                        q_value = regression_info[mod_key][Q_VALUE]

                    if pprpa[phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN].get(METRICS):
                        # TODO - why potentially no metrics? Not generated due to lack of data?
                        pprpa[phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = q_value

                    fisher_output = DEFAULT_FISHER_STATS

                    if mod_key in time_course_fisher_dict:
                        fisher_output = time_course_fisher_dict[mod_key]

                    pprpa[phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS]["Fisher_G"] = fisher_output

            rr.save()

    # calculateProteinOscillationAbundances
    def _calculate_protein_oscillation(self, protein, phospho, replicates, sample_stages, zero_or_log2, statistic_type_name):
        # TODO - gets protein abundances repeatedly, once for each phospho. Inefficient.
        #   Put it outside the loop?
        protein_abundances = self._get_abundances(zero_or_log2, protein=protein)
        phospho_abundances = self._get_abundances(zero_or_log2, phospho=phospho)

        _, stat = self._clear_and_fetch_stats(statistic_type_name, phospho=phospho)
        
        for replicate in replicates:
            # TODO - inefficient
            # TODO - variables poorly named
            protein_replicate = protein_abundances.filter(replicate=replicate)
            phospho_replicate = phospho_abundances.filter(replicate=replicate)

            if not protein_replicate.exists() or not phospho_replicate.exists():
                continue
    
            for sample_stage in sample_stages:
                try:
                    protein_stage = protein_replicate.get(sample_stage=sample_stage)
                    phospho_stage = phospho_replicate.get(sample_stage=sample_stage)
                except:
                    continue

                oscillation = phospho_stage.reading - protein_stage.reading

                Abundance.objects.create(
                    statistic = stat,
                    replicate = replicate,
                    sample_stage = sample_stage,
                    reading = oscillation,
                )


    def _calculate_abundances_metrics(
        self,
        replicates,
        sample_stages,
        protein = None,
        phospho = None,
        with_bugs = False
    ):
        if not protein and not phospho:
            logger.error("Either a protein or a phospho must be provided.")

        # Clear all abundances other than non-average raw
        if protein:
            # Abundance.objects.exclude(
            #     statistic__statistic_type__name=ABUNDANCES_RAW,
            # ).filter(
            #     statistic__protein=protein
            # ).delete()

            Abundance.objects.filter(
                statistic__statistic_type__name=ABUNDANCES_RAW,
                statistic__protein=protein,
                replicate__mean = True
            ).delete()
        else:
            # Abundance.objects.exclude(
            #     statistic__statistic_type__name=ABUNDANCES_RAW
            # ).filter(
            #     statistic__phospho=phospho
            # ).delete()

            Abundance.objects.filter(
                statistic__statistic_type__name=ABUNDANCES_RAW,
                statistic__phospho=phospho,
                replicate__mean = True
            ).delete()

        # firstLevelNormalisationProteomics
        # firstLevelNormalisationPhospho
        # TODO - remove all returned values, e.g. normalised_medians
        normalised_medians = self._calculate_normalised_medians(
            protein, phospho
        )

        # calclog2PalboNormalisation
        arrest_readings = self._calculate_arrest_log2_normalisation(
            protein, phospho
        )

        # calclog2RelativeAbundance
        log2_readings = self._calculate_relative_log2_normalisation(
            protein, phospho
        )

        # normaliseData
        min_max_readings = self._calculate_zero_or_min_normalisation(
            replicates, protein, phospho, False
        )

        # normaliseData
        zero_max_readings = self._calculate_zero_or_min_normalisation(
            replicates, protein, phospho, True
        )

        imputed_readings = self._impute(
            replicates, protein, phospho, 
        )

        raw_averages = (
            self._calculate_means(
                ABUNDANCES_RAW,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs
            )
        )

        normalised_averages = (
            self._calculate_means(
                ABUNDANCES_NORMALISED_MEDIAN,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs
            )
        )

        min_max_averages = (
            self._calculate_means(
                ABUNDANCES_NORMALISED_MIN_MAX,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs
            )
        )

        zero_max_averages = (
            self._calculate_means(
                ABUNDANCES_NORMALISED_ZERO_MAX,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs
            )
        )

        log2_averages = (
            self._calculate_means(
                ABUNDANCES_NORMALISED_LOG2_MEAN,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs
            )
        )

        arrest_averages = (
            self._calculate_means(
                ABUNDANCES_NORMALISED_LOG2_ARREST,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs
            )
        )

        imputed_averages = (
            self._calculate_means(
                ABUNDANCES_IMPUTED,
                protein = protein,
                phospho = phospho,
                with_bugs = with_bugs,
                imputed = True,
            )
        )

        log2_mean_metrics = self._calculate_metrics(
            ABUNDANCES_NORMALISED_LOG2_MEAN,
            replicates,
            sample_stages,
            protein,
            phospho,
        )

        zero_max_mean_metrics = self._calculate_metrics(
            ABUNDANCES_NORMALISED_ZERO_MAX,
            replicates,
            sample_stages,
            protein,
            phospho,
        )

        anova = self._calculate_ANOVA(
            ABUNDANCES_NORMALISED_LOG2_MEAN,
            sample_stages,
            protein,
            phospho
        )

        # # TODO - this is needed in batch processing
        # # relative_log2_readings_by_protein[
        # #     protein
        # # ] = log2_readings



    # TODO - tested
    def _generate_phospho_regression_metrics(self, run, replicates, sample_stages, with_bugs):
        run_results = self._fetch_run_results(run)

        # TODO - this is the reason the fields are in the wrong order. They're not
        #   in the same order when returned from the DB I suspect.
        stages = [s.name for s in sample_stages]

        for rr in run_results:
            # TODO - these names are inconsistent with elsewhere
            phosphoryl_abundances = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]
            protein_abundances = rr.combined_result[PROTEIN_ABUNDANCES]

            if not len(phosphoryl_abundances) or not len(protein_abundances[RAW]):
                continue

            try:
                for phosphosite in phosphoryl_abundances:
                    # TODO - not in the original
                    pappan = phosphoryl_abundances[
                        phosphosite][POSITION_ABUNDANCES][NORMALISED]

                    rep_proteos = []
                    rep_phosphos = []

                    for replicate in replicates:
                        if pappan[LOG2_MEAN].get(replicate.name) is None:
                            continue

                        rep_proteo = {}
                        rep_phospho = {}

                        unordered_rep_proteo = protein_abundances[NORMALISED][LOG2_MEAN][replicate.name]
                        unordered_rep_phospho = pappan[LOG2_MEAN][replicate.name]

                        if len(unordered_rep_proteo) != len(sample_stages) or len(unordered_rep_phospho) != len(sample_stages):
                            continue

                        # Reorder the stages, they come out of the DB in the wrong order
                        for stage in sample_stages:
                            rep_proteo[stage.name] = unordered_rep_proteo[stage.name]
                            rep_phospho[stage.name] = unordered_rep_phospho[stage.name]

                        rep_phosphos.append(rep_phospho)
                        rep_proteos.append(rep_proteo)

                    if len(rep_phosphos) != len(replicates) or len(rep_proteos) != len(replicates):
                        # One of the replicates didn't have all the required data
                        continue

                    # Add the Regressed Phospho Normalised Abundances
                    phosphoryl_abundances[
                        phosphosite][PHOSPHO_REGRESSION] = {ZERO_MAX:{}, LOG2_MEAN:{}}

                    # converting dictionary values to list for both replicates
                    # TODO - study this
                    # TODO - Can this function can be simplified?
                    phospho_Y_list = [val for r_phos in rep_phosphos for val in r_phos.values()]
                    protein_X_list = [val for r_prot in rep_proteos for val in r_prot.values()]

                    # converting list to array
                    prot_x = np.asarray(protein_X_list).reshape(len(protein_X_list), 1)
                    phospho_y = np.asarray(phospho_Y_list).reshape(len(phospho_Y_list), 1)

                    # TODO - a try within a try?
                    try:
                        if None in prot_x or None in phospho_y:
                            continue

                        model = linear_model.LinearRegression().fit(prot_x, phospho_y)
                    except Exception as e:
                        logger.error("_generate_phospho_regression_metrics error")
                        logger.error(prot_x)
                        logger.error(phospho_y)
                        logger.error(e)

                    y_pred = model.predict(prot_x)

                    residuals = (phospho_y - y_pred)

                    # TODO - change this variable name
                    res_dic = {}

                    num_stages = len(sample_stages)

                    # replicates is a list with len len(replicates) * len(sample_stages)
                    #   Chop it up in to bits for each replicate.
                    for i, replicate in enumerate(replicates):
                        res_dic[replicate.name] = {}

                        for j, value in enumerate(residuals[(i*num_stages):((i+1)*num_stages)]):
                            key = stages[j]
                            res_dic[replicate.name][key] = value[0]

                    phospho_regression = phosphoryl_abundances[phosphosite][
                        PHOSPHO_REGRESSION]
                    
                    phospho_regression[LOG2_MEAN] = res_dic

                    phospho_regression[LOG2_MEAN][ABUNDANCE_AVERAGE] = self._calculate_means(
                        phospho_regression[LOG2_MEAN],
                        imputed=False,
                        with_bugs=with_bugs
                    )

                    # Calculate the protein - phospho vector correlation
                    phospho_regression[ZERO_MAX][METRICS] = {}

                    # [0] to get the correlation coefficient, [1] = p-value
                    phospho_regression[ZERO_MAX][METRICS][PROTEIN_PHOSPHO_CORRELATION] = stats.pearsonr(protein_X_list, phospho_Y_list)[0]

                    # curve fold change phosphorylation/curve fold change protein for 0-max
                    phospho_regression[ZERO_MAX][METRICS][PHOSPHO_PROTEIN_CFC_RATIO] = rr.combined_result[PHOSPHORYLATION_ABUNDANCES][
                        phosphosite][METRICS][ZERO_MAX][CURVE_FOLD_CHANGE] / rr.combined_result[METRICS][ZERO_MAX][CURVE_FOLD_CHANGE]

                    # ANOVA
                    phospho_regression[LOG2_MEAN][METRICS] = {}
                    phospho_regression[LOG2_MEAN][METRICS][ANOVA] = self._calculate_ANOVA(phospho_regression[LOG2_MEAN], replicates, sample_stages)
            except Exception as e:
                logger.error("generate_phospho_regression_metrics error")
                logger.error(e)
                logger.error()

            rr.save()

    # TODO - tested
    def _generate_protein_oscillation_metrics(self, project, replicates, sample_stages, with_bugs):
        # run_results = self._fetch_run_results(run)

        # proteins = Protein.objects.filter(project=project).iterator(chunk_size=100)

        # # TODO - could maybe just loop through phosphos here:
        # #   phosphos = Phospho.objects.filter(protein__project=project)
        # #   as checks for values are done in calculate_protein_oscillation
        # # for rr in run_results:
        # for protein in proteins:
        #     # prpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

        #     # Not all samples have both protein and phosphorylation abundances
        #     stat_protein_raw = self._get_statistic(ABUNDANCES_RAW, protein=protein)

        #     print(stat_protein_raw)

        #     # TODO - probably inefficient
        #     if not self._get_abundances(stat_protein_raw).count():
        #         continue

        #     # if not len(prpa) or not len(rr.combined_result[PROTEIN_ABUNDANCES][RAW]):
        #     #     continue

        #     phosphos = Phospho.objects.get(protein=protein)

        phosphos = Phospho.objects.filter(
            protein__project = project
        ).iterator(chunk_size=100)

        # phosphos = Phospho.objects.filter(
        #     protein__project = project
        # )[:100]

        num_phosphos = 0

            # for mod in prpa:
        for phospho in phosphos:
            if not num_phosphos % 1000:
                logger.info(f"Processing oscillation {num_phosphos} {phospho.protein.accession_number} {phospho.mod}")

            num_phosphos += 1
        
            # # Not all samples have both protein and phosphorylation abundances
            # stat_phospho_raw = self._get_statistic(ABUNDANCES_RAW, phospho=phospho)
            # print("+++++ SPR")
            # print(stat_phospho_raw)

            # # TODO - probably inefficient
            # if not self._get_abundances(stat_phospho_raw).count():
            #     continue

            # prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES] = {}
            # prpampoa = prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES]

            for zero_or_log2, statistic_type_name in [
                [ABUNDANCES_NORMALISED_ZERO_MAX, PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX], 
                [ABUNDANCES_NORMALISED_LOG2_MEAN, PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN]
            ]:
                # TODO - should this be passing phosphosite, not mod?
                #   In the original the two seem to be interchangeable to an extent.

                # phospho_oscillations = self._calculat_protein_oscillation(
                #     pr, norm_method, mod, replicates, sample_stages
                # )

                self._calculate_protein_oscillation(
                    phospho.protein, phospho, replicates, sample_stages, zero_or_log2, statistic_type_name
                )

                # prpampoa[norm_method] = phospho_oscillations

                # prpampoa[norm_method][ABUNDANCE_AVERAGE] = self._calculate_means(
                #     prpampoa[norm_method], imputed = False, with_bugs = with_bugs
                # )

                self._calculate_means(
                    statistic_type_name,
                    protein = None,
                    phospho = phospho,
                    with_bugs = with_bugs
                )

                self._calculate_metrics(
                    statistic_type_name,
                    replicates,
                    sample_stages,
                    None,
                    phospho
                )

                    # for norm_method in [ZERO_MAX, LOG2_MEAN]:
                    #     prpampoa[norm_method][METRICS] = self._calculate_metrics(
                    #         prpampoa[norm_method],
                    #         prpampoa[norm_method][ABUNDANCE_AVERAGE],
                    #         replicates,
                    #         sample_stages
                    #     )

                self._calculate_ANOVA(
                    statistic_type_name,
                    sample_stages,
                    None,
                    phospho
                )

                    # anovas = self._calculate_ANOVA(prpampoa[LOG2_MEAN], replicates, sample_stages)

            #         prpampoa[LOG2_MEAN][METRICS][ANOVA] = anovas

            # rr.save()


    # TODO - test
    def _add_q_value(self, info):
        info_df = pd.DataFrame(info)
        info_df = info_df.T

        info_df[Q_VALUE] = stats.false_discovery_control(info_df[P_VALUE])
        # info_df[Q_VALUE] = self.p_adjust_bh(info_df[P_VALUE])

        return info_df.to_dict('index')

    # TODO - lifted from ICR, rewrite
    # TODO - does it really need to be a dataframe? Write tests and change if possible.
    def p_adjust_bh(self, p):
        """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
        p = np.asarray(p, dtype=float)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / np.arange(len(p), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))

        return q[by_orig]


    # getConsensusKinasePred
    def _generate_kinase_predictions(self, run):
        logger.info("Generating kinase predictions")

        run_results = self._fetch_run_results(run)

        for rr in run_results:
            pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            for mod in pprpa:
                # only for certain phosphorylation sites
                # TODO - why? What does the lack of dash mean?

                if mod.find("-") == -1:
                    # Add PepTools Phospho annotations
                    mod_result = pprpa[mod]
                    accession_number = rr.protein.accession_number

                    if PHOSPHO in index_protein_names[accession_number] and mod in index_protein_names[accession_number][PHOSPHO]:
                        mod_result[PEPTOOLS_ANNOTATIONS] = index_protein_names[accession_number][PHOSPHO][mod]

                    mod_result[KINASE_PREDICTION] = {}

                    phospho_kinases_class = self._get_consensus_kinase_pred(accession_number, pprpa[mod][PHOSPHORYLATION_SITE])

                    mod_result[KINASE_PREDICTION][PEPTIDE_SEQ] = phospho_kinases_class[PEPTIDE_SEQ]
                    mod_result[KINASE_PREDICTION][CONSENSUS_MOTIF_MATCH] = phospho_kinases_class[KINASE_MOTIF_MATCH]

            rr.save()

    def _calculate_medians(
        self,
        project,
        replicates,
        sample_stages,
        is_protein,
        with_bugs
    ):
        if is_protein:
            logger.info("Calculating protein medians")
        else:
            logger.info("Calculating phospho medians")

        if is_protein:
            _, stat_prot_med = self._clear_and_fetch_stats(PROTEIN_MEDIAN, project=project)

            abundances = Abundance.objects.filter(
                statistic__protein__project=project,
                statistic__statistic_type__name=ABUNDANCES_RAW
            ).iterator(chunk_size=100)
        else:
            _, stat_prot_med = self._clear_and_fetch_stats(PHOSPHO_MEDIAN, project=project)

            abundances = Abundance.objects.filter(
                statistic__phospho__protein__project=project,
                statistic__statistic_type__name=ABUNDANCES_RAW
            ).iterator(chunk_size=100)

        rep_stage_abundances = {}

        for i, abundance in enumerate(abundances):
            if not i % 10000:
                logger.info(f"Processing protein abundance for median {i}")

            if rep_stage_abundances.get(abundance.replicate) is None:
                rep_stage_abundances[abundance.replicate] = {}

            if rep_stage_abundances[abundance.replicate].get(abundance.sample_stage) is None:
                rep_stage_abundances[abundance.replicate][abundance.sample_stage]  = []

            rep_stage_abundances[abundance.replicate][abundance.sample_stage].append(abundance.reading)

        if is_protein and with_bugs:
            replicate1 = Replicate.objects.get(project=project, name=ICR_ABUNDANCE_REP_1)
            replicate2 = Replicate.objects.get(project=project, name=ICR_ABUNDANCE_REP_2)

            rep_stage_abundances[replicate1] = rep_stage_abundances[replicate2]

        for replicate in replicates:
            if rep_stage_abundances.get(replicate) is None:
                continue

            for sample_stage in sample_stages:
                if rep_stage_abundances[replicate].get(sample_stage) is None:
                    continue

                if not len(rep_stage_abundances[replicate][sample_stage]):
                    logger.error(f"Median with no values (??) {replicate.name} {sample_stage.name}")
                    continue

                median = statistics.median(rep_stage_abundances[replicate][sample_stage])

                Abundance.objects.create(
                    statistic=stat_prot_med,
                    replicate=replicate,
                    sample_stage=sample_stage,
                    reading=median
                )



    def calculate_fisher_g(self, project, replicates, sample_stages, phospho = False, phospho_ab = False, phospho_reg = False):
        logger.info("Calculate Fisher G Statistics")

        time_course_fisher = self._create_abundance_dataframe(project, replicates, sample_stages, phospho, phospho_ab, phospho_reg)
        time_course_fisher = time_course_fisher.dropna()

        # ptest = importr("ptest")

        # for index, row in time_course_fisher.iterrows():
        #     row_z = [i for i in row.tolist()]
        #     z = FloatVector(row_z)

        #     ptestg_res = ptest.ptestg(z, method="Fisher")

        #     g_stat = ptestg_res.rx2("obsStat")[0]
        #     p_value = ptestg_res.rx2("pvalue")[0]
        #     freq = ptestg_res.rx2("freq")

        #     time_course_fisher.loc[index, G_STATISTIC] = g_stat
        #     time_course_fisher.loc[index, P_VALUE] = p_value
        #     time_course_fisher.loc[index, FREQUENCY] = list(freq)

        for index, row in time_course_fisher.iterrows():
            row_z = [i for i in row.tolist()]

            result = ptestg(row_z)

            time_course_fisher.loc[index, G_STATISTIC] = result["obsStat"]
            time_course_fisher.loc[index, P_VALUE] = result["pvalue"]
            time_course_fisher.loc[index, FREQUENCY] = result["freq"]

        q_value = stats.false_discovery_control(time_course_fisher[P_VALUE])

        time_course_fisher[Q_VALUE] = q_value
    
        cols = time_course_fisher.columns
        fisher_cols = [G_STATISTIC, P_VALUE, FREQUENCY, Q_VALUE]
        ab_col = [x for x in cols if x not in fisher_cols]
        time_course_fisher = time_course_fisher.drop(columns=ab_col)
        time_course_fisher_dict = time_course_fisher.to_dict('index')

        return time_course_fisher_dict

    def _add_protein_annotations(self, run):
        logger.info("Adding protein annotations")

        run_results = self._fetch_run_results(run)

        num_proteins = 0

        for rr in run_results:
            if not num_proteins % 1000:
                print(f"Adding annotation {num_proteins} {rr.protein.accession_number}")

            num_proteins += 1

            pan = rr.protein.accession_number

            if index_protein_names.get(pan) is None:
                # TODO - fetch protein info remotely
                # basic_localisation, localisation_keyword = getProteinLocalisation(pan)
                continue

            ipnan = index_protein_names[pan]

            # Not all entries have a value for basic_localisation
            basic_localisation = index_protein_names[pan].get("basic_localisation")
            localisation_keyword = index_protein_names[pan].get("localisation_keyword")

            protein_info = {}

            protein_info["localisation_info"] = {
                "basic_localisation": basic_localisation,
                "localisation_keyword": localisation_keyword
            }

            if "halflife_mean" in ipnan:
                for field in PROTEIN_INFO_FIELDS:
                    protein_info[field] = ipnan[field]     

            rr.combined_result[PROTEIN_INFO] = protein_info

            rr.combined_result[GENE_NAME] = ipnan[GENE_NAME]
            rr.combined_result[PROTEIN_NAME] = ipnan[PROTEIN_NAME]

            # TODO - put in later
            # else:
            #     gene_name, protein_name = getProteinInfo(pan)

            # TODO - consider adding bulk updating
            rr.save()

    # TODO - all records are cleared at the top of _calcualte_abundance_metrics, this is no
    #   longer needed
    def _clear_and_fetch_stats(self, statistic_type_name, project = None, protein = None, phospho = None):
        statistic_type = self._get_statistic_type(statistic_type_name)

        if project:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, project=project)
        elif protein:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, protein=protein)
        elif phospho:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, phospho=phospho)
        else:
            raise Exception(f"_clear_and_fetch_stats needs a project, protein or phospho")

        self._delete_abundances(stat)

        return statistic_type, stat


    # TODO - not needed
    def _fetch_run_results(self, run):
        # Only get the necessary fields to save on memory usage
        return RunResult.objects.only("run", "protein", "combined_result").filter(run=run).iterator(chunk_size=100)

    def _round(self, value):
        return round(value, 4)

    def _get_site_key(self, statistic: Statistic):
        return f"{statistic.phospho.protein.accession_number}_{statistic.phospho.phosphosite}"
    
import numpy as np
from scipy.signal import periodogram
from scipy.special import comb

def fisher_g_test(z):
    """
    Perform Fisher's g-test on a time series `z`.
    Returns:
        g_stat: the observed test statistic
        p_value: p-value for the test
        freq: frequency with max periodogram value
    """
    z = np.asarray(z)
    n = len(z)
    m = (n - 2) // 2 if n % 2 == 0 else (n - 1) // 2

    # Compute the periodogram (uses normalized frequencies)
    freqs, pgram_vals = periodogram(z, scaling='spectrum')

    # Remove Nyquist frequency for even-length input
    if n % 2 == 0:
        pgram_vals = pgram_vals[:-1]

    max_index = np.argmax(pgram_vals)
    g_stat = pgram_vals[max_index] / np.sum(pgram_vals)

    # Compute p-value
    p = int(np.floor(1 / g_stat))
    i_vals = np.arange(1, p + 1)
    terms = comb(m, i_vals) * (-1) ** (i_vals - 1) * (1 - i_vals * g_stat) ** (m - 1)
    p_value = np.sum(terms)

    return np.array([g_stat, p_value, freqs[max_index]])

def ptestg(z):
    """
    Perform periodicity test on input data `z` using Fisher's g-test.
    """
    z = np.atleast_2d(z).T  # Ensure z is a column vector
    n, m = z.shape

    obs_stat = np.zeros(m)
    p_values = np.zeros(m)
    freqs = np.zeros(m)

    for i in range(m):
        result = fisher_g_test(z[:, i])
        obs_stat[i] = result[0]
        p_values[i] = result[1]
        freqs[i] = result[2]

    return {
        "obsStat": obs_stat,
        "pvalue": p_values,
        "freq": freqs,
        "class": "Htest"
    }
