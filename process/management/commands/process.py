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
    PROTEIN_ABUNDANCES_RAW,
    PHOSPHO_READINGS,
    PROTEIN_MEDIAN,

    PROTEIN_ABUNDANCES_RAW,
    PROTEIN_ABUNDANCES_IMPUTED,

    PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX,
    PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
    PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX,
    PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
    PROTEIN_ABUNDANCES_NORMALISED_LOG2_ARREST,

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
            self._calculate_protein_medians(project, replicates, sample_stages, with_bugs)

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
                    protein,
                    with_bugs
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
                        protein,
                        with_bugs
                    )

        # if calculate_phospho_medians or calculate_all:
        #     phospho_medians = self._calculate_phospho_medians(phospho_readings, run)

        # if calculate_phosphos or calculate_all:
        #     self._phospho(project, replicates, phospho_readings, phospho_medians, column_names, phosphos, sample_stages, run, proteins, with_bugs)

        # # Largely batch processing from now on
        # # TODO - some of the parts of the batch processing, e.g. protein oscillation
        # #   metrics, aren't actually batch. Move them out.
        # if calculate_batch or calculate_all:
        #     # TODO - not a batch call, shouldn't be here.
        #     self._add_protein_annotations(run)

        #     self._calculate_phosphorylation_abundances_q_values(run, replicates, sample_stages)

        #     # TODO - figure out how, or if at all, to get this working for SL
        #     #   Actually it no longer works for ICR either, related to
        #     #   P04264 somehow
        #     # self._generate_kinase_predictions(run)

        #     self._calculate_batch_q_value_fisher(run, replicates, sample_stages)

        #     self._add_protein_oscillations(run, replicates, sample_stages, with_bugs)

        #     self._add_phospho_regression(run, replicates, sample_stages, with_bugs)



    def _phospho(self, project, replicates, phospho_readings, phospho_medians, column_names, phosphos, sample_stages, run, proteins, with_bugs: bool):
        logger.info("Processing phosphoproteome")

        # Remove any earlier phospho_result values for this project
        RunResult.objects.filter(run=run).update(phospho_result=None)

        num_processed = 0

        for pr in proteins:
            result = {}
            num_mods = 0

            phospho_readings_filtered = phospho_readings.filter(phospho__protein = pr)

            if not phospho_readings_filtered:
                self._update_phospho_result(run, pr, {})

                continue

            readings = self._format_phospho_readings(phospho_readings_filtered)

            for mod, readings in readings.items():
                num_mods += 1

                self._count_logger(
                    num_mods,
                    10,
                    f"Processing phosphos for {pr.accession_number}, {num_mods}",
                )

                # TODO - is there a more efficient way of doing this?
                #   Maybe prepopulate a dict of phosphos
                phospho = Phospho.objects.get(protein = pr, mod = mod)
                phosphosite = phospho.phosphosite

                result[mod] = {
                    PHOSPHORYLATION_SITE: phosphosite,
                    POSITION_ABUNDANCES: {
                        RAW: {},
                        NORMALISED: {},
                        IMPUTED: {}
                    },
                    METRICS: {}
                }

                self._calculate_abundances_metrics(
                    result[mod],
                    project,
                    replicates,
                    readings,
                    column_names,
                    sample_stages,
                    phospho_medians,
                    POSITION_ABUNDANCES,
                    with_bugs
                )

            run_result, _ = RunResult.objects.get_or_create(
                run=run,
                protein=pr
            )

            run_result.phospho_result = result

            run_result.save()

            num_processed += 1

        logger.info(f"Number of proteins with phospho readings: {num_processed}")

        return result


    def _update_phospho_result(self, run, protein, result):
        run_result, _ = RunResult.objects.get_or_create(
            run=run,
            protein=protein
        )

        run_result.phospho_result = result

        run_result.save()
    

    # TODO - similar to other functionality, consolidate
    def _calculate_phosphorylation_abundances_q_values(self, run, replicates, sample_stages):
        logger.info("Calculating q values for phosphorylation_abundances.")

        # TODO - blank all q_values in DB?

        run_results = self._fetch_run_results(run)

        # TODO - rename this
        prot_anova_info: dict = {}

        for rr in run_results:
            for mod in rr.combined_result[PHOSPHORYLATION_ABUNDANCES]:
                pprpam = rr.combined_result[PHOSPHORYLATION_ABUNDANCES][mod]

                phospho_key = f"{rr.protein.accession_number}_{pprpam['phosphorylation_site']}"

                prot_anova_info[phospho_key] = {
                    P_VALUE: pprpam[METRICS][LOG2_MEAN][ANOVA][P_VALUE]
                }

        # prot_anova_info_df = pd.DataFrame(prot_anova_info).T
        # prot_anova_info_df[Q_VALUE] = self.p_adjust_bh(prot_anova_info_df[P_VALUE])

        # prot_anova_info = prot_anova_info_df.to_dict("index")

        prot_anova_info = self._add_q_value(prot_anova_info)

        run_results = self._fetch_run_results(run)

        for rr in run_results:
            for mod in rr.combined_result[PHOSPHORYLATION_ABUNDANCES]:
                pprpam = rr.combined_result[PHOSPHORYLATION_ABUNDANCES][mod]

                phospho_key = f"{rr.protein.accession_number}_{pprpam['phosphorylation_site']}"

                q_value = 1

                if prot_anova_info.get(phospho_key):
                    q_value = prot_anova_info[phospho_key][Q_VALUE]

                pprpam[METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = q_value

            rr.save()

    def _calculate_batch_q_value_fisher(self, run, replicates, sample_stages):
        logger.info("Calculating q and fisher G values for proteins.")

        # TODO - blank all q_values in DB?

        # Only get the necessary fields to save on memory usage
        run_results_with_log2_mean = RunResult.objects.only("run", "protein", "combined_result").filter(run=run, combined_result__metrics__has_key = "log2_mean")

        fisher_stats = self.calcFisherG(run, replicates, sample_stages)

        # TODO - rename this
        prot_anova_info: dict = {}

        for rr in run_results_with_log2_mean:
            prot_anova_info[rr.protein] = {}
            prot_anova_info[rr.protein][P_VALUE] = rr.combined_result[METRICS][LOG2_MEAN][ANOVA][P_VALUE]

        # prot_anova_info_df = pd.DataFrame(prot_anova_info).T
        # prot_anova_info_df[Q_VALUE] = self.p_adjust_bh(prot_anova_info_df[P_VALUE])

        # prot_anova_info = prot_anova_info_df.to_dict("index")

        prot_anova_info = self._add_q_value(prot_anova_info)

        for rr in run_results_with_log2_mean:
            rr.combined_result[METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = prot_anova_info[rr.protein][Q_VALUE]

            fisher_output = DEFAULT_FISHER_STATS

            if rr.protein.accession_number in fisher_stats:
                fisher_output = fisher_stats[rr.protein.accession_number]

            rr.combined_result[METRICS][LOG2_MEAN][FISHER_G] = fisher_output

            rr.save()



    # TODO - lifted from ICR, rewrite
    # TODO - try to get rid of dataframes if possible
    # TODO - needs a test if to be used, is only used by fisher
    # createAbundanceDf
    def _create_results_dataframe(self, run, replicates, sample_stages, phospho, phospho_ab, phospho_reg):
        # TODO - raw is always False, ignore all raw code
        # TODO - norm_method is always log2_mean, ignore all norm_method code

        run_results = self._fetch_run_results(run)

        abundance_table = {}

        for rr in run_results:
            pan = rr.protein.accession_number

            if phospho:
                for mod in rr.combined_result[PHOSPHORYLATION_ABUNDANCES]:
                    ppam = rr.combined_result[PHOSPHORYLATION_ABUNDANCES][mod]

                    mod_key = f"{pan}_{ppam[PHOSPHORYLATION_SITE]}"

                    if not len(ppam[POSITION_ABUNDANCES]):
                        continue

                    protein_abundances = ppam[POSITION_ABUNDANCES][NORMALISED][LOG2_MEAN]

                    if phospho_ab:
                        if PROTEIN_OSCILLATION_ABUNDANCES not in ppam:
                            continue

                        protein_abundances = ppam[PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN]
                    elif phospho_reg:
                        if PHOSPHO_REGRESSION not in ppam:
                            continue

                        protein_abundances = ppam[PHOSPHO_REGRESSION][LOG2_MEAN]

                    self._build_phospho_abundance_table(abundance_table, replicates, sample_stages, protein_abundances, mod_key)
            else:
                # Not all proteins have protein results, some are phospho only
                if not rr.combined_result[PROTEIN_ABUNDANCES][NORMALISED].get(LOG2_MEAN):
                    continue

                pprpanlm = rr.combined_result[PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN]

                abundance_table[pan] = {}

                for replicate in replicates:
                    for stage in sample_stages:
                        if not pprpanlm[replicate.name].get(stage.name):
                            continue

                        rep_stage_name = f"{replicate.name}_{stage.name}"

                        abundance_table[pan][rep_stage_name] = pprpanlm[
                            replicate.name
                        ][stage.name]

        time_course_abundance_df = pd.DataFrame(abundance_table)
        time_course_abundance_df = time_course_abundance_df.T

        new_cols = []

        for sample_stage in sample_stages:
            for replicate in replicates:
                new_cols.append(f"{replicate.name}_{sample_stage.name}")

        # Rearrange the column order so replicates are near eatch other
        try:
            time_course_abundance_df = time_course_abundance_df[new_cols]
        except Exception as e:
            logger.error("DF FAILED")
            logger.error(phospho)
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
    def _calculate_ANOVA(self, statistic_type_name, protein, replicates, sample_stages):
        # Defaults if not enough replicates
        p_value = 1
        f_statistic = 1

        stat_log2_mean = self._get_statistic(
            statistic_type_name,
            protein=protein
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
            logger.error("Error in _calculate_anova")
            logger.error(e)

    # TODO - lifted from ICR, change names
    def _calculate_metrics(
        self,
        statistic_type_name: str,
        protein: Protein,
        replicates,
        sample_stages
    ):
        metrics = {}

        stat = self._get_statistic(statistic_type_name, protein=protein)
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

            if not readings.get(abfr.replicate.name):
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

            # if we have info for the protein in at least 2 replicates
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

    # TODO - what is this for? Why does it work on PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX?
    def _impute(
        # TODO - all these pr types are wrong, and also probably bad variable names
        self,
        protein: Protein,
        replicates: QuerySet[Replicate],
    ):
        _, stat_imputed = self._clear_and_fetch_stats(
            PROTEIN_ABUNDANCES_IMPUTED,
            protein = protein
        )

        abundances = self._get_abundances(
            PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX,
            protein = protein
        ).order_by(
            "sample_stage__rank"
        )

        abundances_by_rep = {}

        for abundance in abundances:
            if not abundances_by_rep.get(abundance.replicate):
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
                        # TODO - find out why this calculation
                        # TODO - name variables better
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

        # replicates_by_name: dict = {}
        # column_names_by_replicate: dict = {}

        # # TODO - is this needed now we no longer use Replicate objects as keys?
        # for replicate in replicates:
        #     replicates_by_name[replicate.name] = replicate
        #     column_names_by_replicate[replicate.name] = []

        # # TODO - and this?
        # for column_name in column_names:
        #     column_names_by_replicate[column_name.replicate.name].append(
        #         column_name.sample_stage.name
        #     )

        # imputed_readings: dict = {}

        # for replicate_name in readings.keys():
        #     imputed_readings[replicate_name] = {}

        #     abundances_dict = readings[replicate_name]
        #     abundances = list(abundances_dict.values())

        #     stage_names = column_names_by_replicate[replicate_name]

        #     for idx, stage_name in enumerate(stage_names):
        #         # Default value, should never be used
        #         value = 0

        #         if abundances_dict.get(stage_name) is not None:
        #             value = abundances_dict[stage_name]
        #         else:
        #             last = None
        #             next = None

        #             # TODO - isn't there a better way to iterate?
        #             for offset in range(1, len(stage_names)):
        #                 prev_idx = idx - offset
        #                 if prev_idx < 0:
        #                     # Gone before the beginning of the list, give up
        #                     break

        #                 prev_stage_name = stage_names[prev_idx]

        #                 if abundances_dict.get(prev_stage_name) is not None:
        #                     last = (offset, abundances_dict[prev_stage_name])
        #                     # last = abundances_dict[prev_stage_name]
        #                     break

        #             for offset in range(1, len(abundances)):
        #                 # Look forward
        #                 # TODO - this seems to loop back to the beginning. Is that right?
        #                 next_idx = (idx + offset) % len(abundances)
        #                 next_stage_name = stage_names[next_idx]

        #                 if abundances_dict.get(stage_name) is not None:
        #                     next = (offset, abundances[next_stage_name])
        #                     # next = abundances[next_stage_name]
        #                     break

        #             if last and next:
        #                 # Linear imputation between nearest timepoints
        #                 # TODO - find out why this calculation
        #                 # TODO - name variables better
        #                 d1, a1 = last
        #                 d2, a2 = next
        #                 step_height = (a1 - a2) / (d1 + d2)
        #                 value = d2 * step_height + a2

        #         # imputed_protein_readings[protein][replicate_name][stage_name] = self._round(float(value))
        #         # TODO - for some reason ICR rounds to 2, not 4. What to do?
        #         imputed_readings[replicate_name][stage_name] = round(float(value), 2)

        # return imputed_readings

    def _calculate_zero_or_min_normalisation(self, protein: Protein, replicates, zero_min=False):
        if zero_min:
            _, stat = self._clear_and_fetch_stats(
                PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX,
                protein = protein
            )

            abundances = self._get_abundances(
                PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
                protein = protein
            )
        else:
            _, stat = self._clear_and_fetch_stats(
                PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX,
                protein = protein
            )

            abundances = self._get_abundances(
                PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
                protein = protein
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

        #         level_two_normalised_readings[replicate_name][stage_name] = self._round(
        #             abundance_normalised
        #         )

        # level_two_normalised_readings: dict = {}

        # for replicate_name in readings:
        #     level_two_normalised_readings[replicate_name] = {}

        #     min_value = 0
        #     max_value = 0

        #     abundances = readings[replicate_name]

        #     abundance_values_non_null = [
        #         val for val in abundances.values() if val is not None
        #     ]

        #     if len(abundance_values_non_null) != 0:
        #         if not zero_min:
        #             min_value = min(abundance_values_non_null)

        #         max_value = max(abundance_values_non_null)

        #     for stage_name, abundance in abundances.items():
        #         denominator = max_value - min_value
        #         if abundance is None or denominator == 0:
        #             abundance_normalised = 0.5
        #         else:
        #             abundance_normalised = (abundance - min_value) / denominator

        #         level_two_normalised_readings[replicate_name][stage_name] = self._round(
        #             abundance_normalised
        #         )

        # return level_two_normalised_readings

    # calclog2RelativeAbundance
    def _calculate_relative_log2_normalisation(self, protein: Protein):
        _, stat_normalised_log2_mean = self._clear_and_fetch_stats(
            PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein = protein
        )

        abundances = self._get_abundances(
            PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
            protein = protein
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

        # for replicate_name in readings:
        #     log2_abundances[replicate_name] = {}

        #     for stage_name in readings[replicate_name]:
        #         log2 = None
        #         reading = readings[replicate_name][stage_name]

        #         if reading is not None:
        #             log2 = math.log2(reading)

        #         log2_abundances[replicate_name][stage_name] = log2

        # total_abundances = 0
        # total_lengths = 0

        # log2_normalised_readings: dict = {}

        # for replicate_name in log2_abundances:
        #     for stage_name in log2_abundances[replicate_name]:
        #         if log2_abundances[replicate_name][stage_name] is not None:
        #             total_abundances += log2_abundances[replicate_name][stage_name]
        #             total_lengths += 1

        #     mean = None

        #     if total_lengths != 0:
        #         mean = total_abundances / total_lengths
        #     # TODO - is mean is None then the loop below can be simplified

        #     log2_normalised_readings = {}

        #     for replicate_name in readings:
        #         log2_normalised_readings[replicate_name] = {}

        #         for stage_name in readings[replicate_name]:
        #             normalised_abundance = None

        #             if log2_abundances[replicate_name].get(stage_name) is not None and mean is not None:
        #                 normalised_abundance = self._round(
        #                     log2_abundances[replicate_name][stage_name] - mean
        #                 )

        #             log2_normalised_readings[replicate_name][
        #                 stage_name
        #             ] = normalised_abundance

        # return log2_normalised_readings

    def _calculate_arrest_log2_normalisation(self, protein: Protein):
        # TODO - is ARRESTING_AGENT the wrong name?
        ARRESTING_AGENT = "Nocodozole"

        # TODO - this is a hack, maybe add the field to the Project model?
        if protein.project.name == "ICR":
            ARRESTING_AGENT = "Palbo"

        _, stat_normalised_log2_arrest = self._clear_and_fetch_stats(
            PROTEIN_ABUNDANCES_NORMALISED_LOG2_ARREST,
            protein = protein
        )

        abundances = self._get_abundances(
            PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
            protein = protein
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


    def _format_phospho_readings_full(self, phospho_readings):
        logger.info(
            "Converting phospho_readings QuerySet into dict by protein, mod, replicate and stage name"
        )

        readings_by_rep_stage: dict = {}

        mod_no = 0

        for phospho_reading in phospho_readings:
            mod_no += 1
            self._count_logger(
                mod_no,
                10000,
                f"Formatting phospho for {mod_no}, protein {phospho_reading.phospho.protein.accession_number} mod {phospho_reading.phospho.mod}",
            )

            # TODO - what about Nones? Will there be any here? Check the import script.
            reading = phospho_reading.reading

            protein = phospho_reading.phospho.protein
            mod = phospho_reading.phospho.mod
            replicate_name = phospho_reading.column_name.replicate.name
            stage_name = phospho_reading.column_name.sample_stage.name

            if not readings_by_rep_stage.get(protein):
                readings_by_rep_stage[protein] = {}

            if not readings_by_rep_stage[protein].get(mod):
                readings_by_rep_stage[protein][mod] = {}

            if not readings_by_rep_stage[protein][mod].get(replicate_name):
                readings_by_rep_stage[protein][mod][replicate_name] = {}

            readings_by_rep_stage[protein][mod][replicate_name][stage_name] = reading

        return readings_by_rep_stage

    def _format_phospho_readings(self, phospho_readings):
        readings_by_rep_stage: dict = {}

        for phospho_reading in phospho_readings:
            # TODO - what about Nones? Will there be any here? Check the import script.
            reading = phospho_reading.reading

            mod = phospho_reading.phospho.mod
            replicate_name = phospho_reading.column_name.replicate.name
            stage_name = phospho_reading.column_name.sample_stage.name

            if not readings_by_rep_stage.get(mod):
                readings_by_rep_stage[mod] = {}

            if not readings_by_rep_stage[mod].get(replicate_name):
                readings_by_rep_stage[mod][replicate_name] = {}

            readings_by_rep_stage[mod][replicate_name][stage_name] = reading

        return readings_by_rep_stage

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
        return self._get_abundances(PROTEIN_ABUNDANCES_RAW, protein=protein)

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

    def _calculate_normalised_medians(self, protein):
        _, stat_normalised_medians = self._clear_and_fetch_stats(
            PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
            protein = protein
        )

        stat_protein_medians = self._get_statistic(PROTEIN_MEDIAN, project = protein.project)

        # TODO - get rid of _get_protein_readings? It may not be that useful
        protein_reading = self._get_protein_readings(protein=protein)

        for prr in protein_reading:
            reading = prr.reading

            # TODO - is this if statement needed?
            if reading is not None:
                median = self._get_abundance(
                    stat_protein_medians, 
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
        protein: Protein,
        with_bugs: bool,
        imputed: bool = False,
    ):
        readings: dict = {}

        abundances = self._get_abundances(statistic_type_name, protein=protein)

        for abundance in abundances:
            sample_stage = abundance.sample_stage

            if not readings.get(sample_stage):
                readings[sample_stage] = []

            if with_bugs and not imputed:
                # We throw away the second reading
                # TODO - how will this behave towards None?
                if len(readings[sample_stage]) == 1:
                    continue

            if abundance.reading is not None:
                readings[sample_stage].append(abundance.reading)

        stat = self._get_statistic(statistic_type_name, protein=protein)

        # Delete all mean abundances for this statistic
        Abundance.objects.filter(statistic=stat, replicate__mean=True).delete()

        mean_replicate = Replicate.objects.get(project=protein.project, mean=True)

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



    def _calculate_phospho_medians(
        self,
        phospho_readings: dict,
        run
    ):
        # TODO - get rid of _format_phospho_readings, or at least put it in the loop
        readings = self._format_phospho_readings_full(phospho_readings)

        phospho_medians = {}

        for protein in readings.keys():
            for mod in readings[protein].keys():
                for replicate_name in readings[protein][mod].keys():
                    if not phospho_medians.get(replicate_name):
                        phospho_medians[replicate_name] = {}

                    for column_name in readings[protein][mod][replicate_name].keys():
                        if not phospho_medians[replicate_name].get(column_name):
                            phospho_medians[replicate_name][column_name] = []

                        reading = readings[protein][mod][replicate_name][column_name]

                        if reading is None:
                            continue

                        phospho_medians[replicate_name][column_name].append(reading)

        for replicate_name in phospho_medians.keys():
            for column_name in phospho_medians[replicate_name].keys():
                median = statistics.median(phospho_medians[replicate_name][column_name])

                phospho_medians[replicate_name][column_name] = median

        run, _ = Run.objects.get_or_create(
            project=run.project, with_bugs=run.with_bugs
        )

        run.phospho_medians = json.dumps(phospho_medians)

        run.save()

        return phospho_medians

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
        #   first removed some of the R code from calcFisherG.
        #   It's something to do with P04264, one of the duplicated proteins
        #   for ICR. It's the reason I made it like get_or_create in import_proteo
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



    # TODO - lifted from ICR
    # addProteinOscillations
    def _add_protein_oscillations(self, run, replicates, sample_stages, with_bugs):
        # TODO - check all logging statements for similarity to ICR
        logger.info("Adding Protein Oscillation Normalised Abundances")

        self._generate_protein_oscillation_metrics(run, replicates, sample_stages, with_bugs)

        # TODO - have this put the values directly into RR?
        time_course_fisher_dict = self.calcFisherG(run, replicates, sample_stages, phospho = True, phospho_ab = True)

        ps_and_qs = {}

        run_results = self._fetch_run_results(run)

        num_proteins = 0

        for rr in run_results:
            pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            self._count_logger(
                num_proteins,
                1000,
                f"Calculating protein oscillation {num_proteins} {rr.protein.accession_number}",
            )

            num_proteins += 1

            for mod in pprpa:
                if pprpa[mod].get(PROTEIN_OSCILLATION_ABUNDANCES):
                    phospho_key = f"{rr.protein.accession_number}_{pprpa[mod][PHOSPHORYLATION_SITE]}"

                    if phospho_key not in ps_and_qs:
                        ps_and_qs[phospho_key] = {}

                    ps_and_qs[phospho_key][P_VALUE] = pprpa[mod][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][P_VALUE]

        # TODO - what's this for?
        ps_and_qs = self._add_q_value(ps_and_qs)

        run_results = self._fetch_run_results(run)

        # TODO - tidy up, two loops not necessary?
        for rr in run_results:
            pprpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            # TODO - why check for raw?
            if not len(pprpa) or not len(rr.combined_result[PROTEIN_ABUNDANCES][RAW]):
                continue

            for mod in pprpa:
                if not pprpa[mod].get(PROTEIN_OSCILLATION_ABUNDANCES):
                    continue

                mod_key = f"{rr.protein.accession_number}_{pprpa[mod][PHOSPHORYLATION_SITE]}"

                q_value = 1

                if mod_key in ps_and_qs:
                    q_value = ps_and_qs[mod_key][Q_VALUE]

                pprpa[mod][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = q_value

                # Fisher
                # TODO - tidy this and others like it
                fisher_output = DEFAULT_FISHER_STATS

                if mod_key in time_course_fisher_dict:
                    fisher_output = time_course_fisher_dict[mod_key]

                pprpa[mod][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS]["Fisher_G"] = fisher_output

            rr.save()



    # TODO - lifted from ICR
    # addPhosphoRegression
    def _add_phospho_regression(self, run, replicates, sample_stages, with_bugs):
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

        self._generate_phospho_regression_metrics(run, replicates, sample_stages, with_bugs)

        # Fisher G Statistic - Phospho
        time_course_fisher_dict = self.calcFisherG(run, replicates, sample_stages, phospho = True, phospho_ab = False, phospho_reg = True)

        regression_info = {}

        run_results = self._fetch_run_results(run)

        num_proteins = 0

        # TODO - looks very close to the code in protein oscillation. Make a function?
        for rr in run_results:
            self._count_logger(
                num_proteins,
                1000,
                f"Adding phospho regression {num_proteins} for {rr.protein.accession_number}",
            )

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
    # TODO - tested
    def _calculate_protein_oscillation(self, result, norm_method, phosphosite, replicates, sample_stages):
        '''
        Iterates through the replicates and stages for
            results_single[PROTEIN_ABUNDANCES][NORMALISED][norm_method]
            result[PHOSPHORYLATION_ABUNDANCES][
                        phosphosite][POSITION_ABUNDANCES][NORMALISED][norm_method]

            subtracts one from the other, puts them in a dict by replicate and
            sample stage and returns it.

            It's basically how much protein and phospho values differ for each normalised
            sample stage.
        '''

        phospho_oscillations = {}

        rspannm = result[PROTEIN_ABUNDANCES][NORMALISED][norm_method]
        rspappann = result[PHOSPHORYLATION_ABUNDANCES][
                        phosphosite][POSITION_ABUNDANCES][NORMALISED][norm_method]

        for replicate in replicates:
            if replicate.name in rspannm and replicate.name in rspappann:
                protein_oscillation_abundances = {}

                protein_normed = rspannm[replicate.name]
                phospho_normed = rspappann[replicate.name]

                try:
                    for sample_stage in sample_stages:
                        if not protein_normed.get(sample_stage.name) or not phospho_normed.get(sample_stage.name):
                            continue

                        protein_oscillation_abundances[sample_stage.name] = (
                            phospho_normed[sample_stage.name]
                            - protein_normed[sample_stage.name]
                        )
                except Exception as e:
                    logger.error("_calculate_protein_oscillation error")
                    logger.error(protein_normed)
                    logger.error(phospho_normed)
                    logger.error(e)
                    exit()

                phospho_oscillations[replicate.name] = protein_oscillation_abundances

        return phospho_oscillations

    def _calculate_abundances_metrics(
        self,
        replicates,
        sample_stages,
        protein,
        with_bugs
    ):
        # Clear any raw average abundances set by a previous run
        Abundance.objects.filter(
            statistic__statistic_type__name=PROTEIN_ABUNDANCES_RAW,
            statistic__protein=protein,
            replicate__mean = True
        ).delete()

        # firstLevelNormalisationProteomics
        # firstLevelNormalisationPhospho
        # TODO - remove all returned values, e.g. normalised_medians
        normalised_medians = self._calculate_normalised_medians(
            protein
        )

        # calclog2PalboNormalisation
        arrest_readings = self._calculate_arrest_log2_normalisation(
            protein
        )

        # calclog2RelativeAbundance
        log2_readings = (
            self._calculate_relative_log2_normalisation(protein)
        )

        # normaliseData
        min_max_readings = self._calculate_zero_or_min_normalisation(
            protein, replicates
        )

        # normaliseData
        zero_max_readings = self._calculate_zero_or_min_normalisation(
            protein, replicates, True
        )

        imputed_readings = self._impute(
            protein, replicates
        )

        raw_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_RAW, protein, with_bugs
            )
        )

        normalised_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_NORMALISED_MEDIAN, protein, with_bugs
            )
        )

        min_max_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX, protein, with_bugs
            )
        )

        zero_max_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX, protein, with_bugs
            )
        )

        log2_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN, protein, with_bugs
            )
        )

        arrest_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_NORMALISED_LOG2_ARREST, protein, with_bugs
            )
        )

        imputed_averages = (
            self._calculate_means(
                PROTEIN_ABUNDANCES_IMPUTED, protein, with_bugs, imputed=True
            )
        )

        log2_mean_metrics = self._calculate_metrics(
            PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein,
            replicates,
            sample_stages,
        )

        zero_max_mean_metrics = self._calculate_metrics(
            PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX,
            protein,
            replicates,
            sample_stages,
        )

        anova = self._calculate_ANOVA(
            PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein,
            replicates,
            sample_stages
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
                        if not pappan[LOG2_MEAN].get(replicate.name):
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
    def _generate_protein_oscillation_metrics(self, run, replicates, sample_stages, with_bugs):
        run_results = self._fetch_run_results(run)

        for rr in run_results:
            prpa = rr.combined_result[PHOSPHORYLATION_ABUNDANCES]

            # Not all samples have both protein and phosphorylation abundances
            if not len(prpa) or not len(rr.combined_result[PROTEIN_ABUNDANCES][RAW]):
                continue

            for mod in prpa:
                prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES] = {}
                prpampoa = prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES]

                for norm_method in [ZERO_MAX, LOG2_MEAN]:
                    # TODO - should this be passing phosphosite, not mod?
                    #   In the original the two seem to be interchangeable to an extent.
                    phospho_oscillations = self._calculate_protein_oscillation(
                        rr.combined_result, norm_method, mod, replicates, sample_stages
                    )

                    prpampoa[norm_method] = phospho_oscillations

                    # TODO - split this out into its own variable? That's done with other average
                    prpampoa[norm_method][ABUNDANCE_AVERAGE] = self._calculate_means(
                        prpampoa[norm_method], imputed = False, with_bugs = with_bugs
                    )

                # TODO - check all comments to make sure they're non-ICR
                # TODO - why only the check for LOG2_MEAN? Why not for both?
                if len(prpampoa[LOG2_MEAN]) > 1:
                    for norm_method in [ZERO_MAX, LOG2_MEAN]:
                        prpampoa[norm_method][METRICS] = self._calculate_metrics(
                            prpampoa[norm_method],
                            prpampoa[norm_method][ABUNDANCE_AVERAGE],
                            replicates,
                            sample_stages
                        )

                    anovas = self._calculate_ANOVA(prpampoa[LOG2_MEAN], replicates, sample_stages)

                    prpampoa[LOG2_MEAN][METRICS][ANOVA] = anovas

            rr.save()


    # TODO - is this worth being a function at all?
    def _build_phospho_abundance_table(self, abundance_table, replicates, sample_stages, protein_abundances, mod_key):
        for replicate in replicates:
            if not abundance_table.get(mod_key):
                abundance_table[mod_key] = {}

            # TODO - not in original, why here?
            if not protein_abundances.get(replicate.name):
                continue

            # TODO - check all these loops have the same variable name
            for sample_stage in sample_stages:
                # rep = "_".join(abundance_rep.split("_", 1)[1].split("_")[:2])
                if not protein_abundances[replicate.name].get(sample_stage.name):
                    continue

                rep_timepoint = f"{replicate.name}_{sample_stage.name}"
                abundance_table[mod_key][rep_timepoint] = protein_abundances[replicate.name][sample_stage.name]

    # TODO - test
    def _add_q_value(self, info):
        info_df = pd.DataFrame(info)
        info_df = info_df.T

        info_df[Q_VALUE] = stats.false_discovery_control(p)

        return info_df.to_dict('index')

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

    def _calculate_protein_medians(self, project, replicates, sample_stages, with_bugs):
        logger.info("Calculating protein medians")

        _, stat_prot_med = self._clear_and_fetch_stats(PROTEIN_MEDIAN, project=project)

        abundances = Abundance.objects.filter(
            statistic__protein__project=project,
            statistic__statistic_type__name=PROTEIN_ABUNDANCES_RAW
        ).iterator(chunk_size=100)

        rep_stage_abundances = {}

        for i, abundance in enumerate(abundances):
            if not i % 10000:
                logger.info(f"Processing abundance for median {i}")

            if not rep_stage_abundances.get(abundance.replicate):
                rep_stage_abundances[abundance.replicate] = {}

            if not rep_stage_abundances[abundance.replicate].get(abundance.sample_stage):
                rep_stage_abundances[abundance.replicate][abundance.sample_stage]  = []

            rep_stage_abundances[abundance.replicate][abundance.sample_stage].append(abundance.reading)

        if with_bugs:
            replicate1 = Replicate.objects.get(project=project, name="abundance_rep_1")
            replicate2 = Replicate.objects.get(project=project, name="abundance_rep_2")

            rep_stage_abundances[replicate1] = rep_stage_abundances[replicate2]

        for replicate in replicates:
            for sample_stage in sample_stages:
                if not rep_stage_abundances[replicate].get(sample_stage):
                    logger.error(f"Median with no sample stage (??) {replicate.name} {sample_stage.name}")

                if not len(rep_stage_abundances[replicate][sample_stage]):
                    logger.error(f"Median with no values (??) {replicate.name} {sample_stage.name}")

                median = statistics.median(rep_stage_abundances[replicate][sample_stage])

                Abundance.objects.create(statistic=stat_prot_med, replicate=replicate, sample_stage=sample_stage, reading=median)


    # TODO - rename
    def calcFisherG(self, run, replicates, sample_stages, phospho = False, phospho_ab = False, phospho_reg = False):
        logger.info("Calculate Fisher G Statistics")

        time_course_fisher = self._create_results_dataframe(run, replicates, sample_stages, phospho, phospho_ab, phospho_reg)
        time_course_fisher = time_course_fisher.dropna()

        ptest = importr("ptest")

        for index, row in time_course_fisher.iterrows():
            row_z = [i for i in row.tolist()]
            z = FloatVector(row_z)

            ptestg_res = ptest.ptestg(z, method="Fisher")

            g_stat = ptestg_res.rx2("obsStat")[0]
            p_value = ptestg_res.rx2("pvalue")[0]
            freq = ptestg_res.rx2("freq")

            time_course_fisher.loc[index, G_STATISTIC] = g_stat
            time_course_fisher.loc[index, 'p_value'] = p_value
            time_course_fisher.loc[index, FREQUENCY] = list(freq)

        q_value = self.p_adjust_bh(time_course_fisher['p_value'])

        time_course_fisher['q_value'] = q_value
    
        cols = time_course_fisher.columns
        fisher_cols = [G_STATISTIC,'p_value', FREQUENCY, 'q_value']
        ab_col = [x for x in cols if x not in fisher_cols]
        time_course_fisher = time_course_fisher.drop(columns=ab_col)
        time_course_fisher_dict = time_course_fisher.to_dict('index')

        return time_course_fisher_dict

    def _add_protein_annotations(self, run):
        logger.info("Adding protein annotations")

        run_results = self._fetch_run_results(run)

        num_proteins = 0

        for rr in run_results:
            self._count_logger(
                num_proteins,
                1000,
                f"Adding annotation {num_proteins} {rr.protein.accession_number}",
            )

            num_proteins += 1

            pan = rr.protein.accession_number

            if not index_protein_names.get(pan):
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

    # TODO - not needed
    def _count_logger(self, i: int, step: int, output: str):
        if i % step == 0:
            logger.info(output)

    def _round(self, value):
        return round(value, 4)
