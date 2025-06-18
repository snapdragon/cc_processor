# TODO - q_value is wrong
#      "ANOVA": {
#        "q_value": 0.6470752162303733
#
#      "ANOVA": {
#        "q value": 0.5788203095177293

# TODO - move this out of _proteo, it is for multiple proteins
#   anovas[protein] = self._calculate_ANOVA(log2_readings)


# TODO - store phospho results in the DB


# TODO - make _proteo iterate by protein and not prebuild the dict
# TODO - put kinase prediction back in. Make them a flag?
# TODO - put PepTools_annotations in
# TODO - phosphorylation metrics log2_mean anova is wrong. Or maybe proteo is. Investigate.
# TODO - ICR doesn't have phospho_regression (??)
# TODO - get rid of 'abundance_average' as an output field name
# TODO - put peptide_abundances in SL?
# TODO - most of the warnings are from _calculate_metrics



import re
import json
import logging
import math
import statistics
import copy
import requests

import numpy as np
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from django.db.models.query import QuerySet
from scipy.signal import periodogram
from scipy.stats import f_oneway, moment
from sklearn.metrics import r2_score
from sklearn import linear_model
from scipy import stats

from process.models import ColumnName, Project, Phospho, Protein, ProteinReading, Replicate, PhosphoReading, SampleStage, Run, RunResult
from process.constants import (PROTEIN_LIMITS,
    FOCUS_PROTEIN_ACCESSION_NUMBER,
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
    TOTAL_PROTEIN_INDEX_FILE)

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
            "--with-bugs",
            help="Run with the original ICR bugs",
            action="store_true",
        )
        parser.add_argument(
            "--limit-proteins",
            help="Limit processing to a subset of proteins",
            action="store_true"
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
            help="Calculate and store protein results",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-phospho-medians",
            help="Calculate and store phospho medians",
            action="store_true"
        )
        parser.add_argument(
            "--calculate-all",
            help="Calculate and store all values",
            action="store_true"
        )

    def handle(self, *args, **options):
        project_name = options["project"]
        with_bugs = options["with_bugs"]
        limit_proteins = options["limit_proteins"]
        calculate_protein_medians = options["calculate_protein_medians"]
        calculate_proteins = options["calculate_proteins"]
        calculate_phospho_medians = options["calculate_phospho_medians"]
        calculate_phosphos = options["calculate_phosphos"]
        calculate_all = options["calculate_all"]

        if not project_name:
            raise Exception("You must provide a project name.")

        if with_bugs and project_name != "ICR":
            raise CommandError("Only an ICR project can run --with-bugs")

        logger.info(f"Processing for project {project_name}, with bugs {with_bugs}")

        # TODO - return if no calculate flags are set

        project = Project.objects.get(name=project_name)
        replicates = Replicate.objects.filter(project=project)
        column_names = ColumnName.objects.filter(replicate__project=project)

        run, _ = Run.objects.get_or_create(project=project, limit_proteins=limit_proteins, with_bugs=with_bugs)

        protein_readings = None
        phospho_readings = None
        phosphos = None
        sample_stages = SampleStage.objects.filter(project=project).order_by('rank')

        if limit_proteins:
            logger.info("Limiting proteins.")
            protein_readings = ProteinReading.objects.filter(
                column_name__replicate__project=project, protein__accession_number__in = PROTEIN_LIMITS
            )
            phospho_readings = PhosphoReading.objects.filter(
                phospho__protein__project=project, phospho__protein__accession_number__in = PROTEIN_LIMITS
            )
            phosphos = Phospho.objects.filter(
                protein__project=project, protein__accession_number__in = PROTEIN_LIMITS
            )
        else:
            protein_readings = ProteinReading.objects.filter(
                column_name__replicate__project=project
            )
            phospho_readings = PhosphoReading.objects.filter(phospho__protein__project=project)
            phosphos = Phospho.objects.filter(protein__project=project)

        FOCUS_PROTEIN = Protein.objects.get(
            accession_number=FOCUS_PROTEIN_ACCESSION_NUMBER, project__name=project.name
        )

        if (calculate_protein_medians or calculate_all) and not limit_proteins:
            protein_medians = self._calculate_protein_medians(run, replicates, protein_readings, column_names, with_bugs)
        else:
            protein_medians = self._fetch_protein_medians(run)

        if calculate_proteins or calculate_all:
            self._proteo(project, replicates, protein_readings, protein_medians, column_names, sample_stages, limit_proteins, run, with_bugs)

        if (calculate_phospho_medians or calculate_all) and not limit_proteins:
            phospho_medians = self._calculate_phospho_medians(phospho_readings, run)
        else:
            phospho_medians = self._fetch_phospho_medians(run)

        if calculate_phosphos or calculate_all:
            self._phospho(project, replicates, phospho_readings, phospho_medians, column_names, phosphos, sample_stages, run, with_bugs)

        # self._merge_phospho_with_proteo(results, phospho_results)

        # self._add_protein_oscillations(results, replicates, sample_stages, with_bugs)

        # results = self._add_phospho_regression(results, replicates, sample_stages, with_bugs)

        # self._save_data(results[FOCUS_PROTEIN], f"{FOCUS_PROTEIN_ACCESSION_NUMBER}.json", False)
        # self._save_data(results, f"{project.name}_limit_protein_{limit_proteins}.json", True)


    def _phospho(self, project, replicates, phospho_readings, phospho_medians, column_names, phosphos, sample_stages, run, with_bugs: bool):
        logger.info("Processing phosphoproteome")

        # Remove any earlier phospho_result values for this project
        RunResult.objects.filter(run=run).update(phospho_result=None)

        # TODO - make this work one by one
        raw_readings = self._format_phospho_readings(phospho_readings)

        # TODO - should this be elsewhere? Or individual queries per loop?
        phosphosites = self._get_phosphosites(phosphos)

        num_proteins = 0

        for protein in raw_readings.keys():
            result = {}

            for mod, readings in raw_readings[protein].items():
                num_proteins += 1

                self._count_logger(
                    num_proteins,
                    10,
                    f"Processing for {num_proteins}, {protein.accession_number}",
                )

                result[mod] = {
                    PHOSPHORYLATION_SITE: phosphosites[mod],
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
                protein=protein
            )

            run_result.phospho_result = result

            run_result.save()

            if protein.accession_number == FOCUS_PROTEIN_ACCESSION_NUMBER:
                self._save_data(result, "Q0966_phospho_only.json", False)

        # TODO - is this batch or single? I'm guessing batch.
        # self._generate_kinase(raw_readings, results)

        return result


    # TODO - rename to protein_medians
    def _proteo(
        self, project, replicates, protein_readings, medians, column_names, sample_stages, limit_proteins, run, with_bugs: bool
    ):
        """
        Does all the required calculations. The steps are:

        TODO - put the descriptions below in the function comments, not here, along with
            data structures.

        1) Calculate the median for each replicate for each stage.
            In other words, get all the values for each column and find the median.

        2) Calculate the mean for each protein for each stage across replicates.
            So for proteins 'ABCD', 'EFGH' with replicates 'One', 'Two' and stages '1h', '2h',
            it will get the abundances for 'ABCD 1h' for each of the replicates, then take
            the mean.

        3) Normalise all abundances. This is done by dividing each abundance by the median
            for its column, then averaging them across replicates.

        TODO - finish this list
        """
        logger.info("Processing proteome")

        # Remove any earlier protein_result values for this project
        RunResult.objects.filter(run=run).update(protein_result=None)

        # N.B. protein_readings_by_rep_stage is not the same structure as protein_readings.
        #   protein_readings is just a list of ProteinReading objects. normalised_protein_readings
        #   is a dict with Protein object keys. The values is a dict of replicate name keys
        #   with a dict of sample stage names and abundances as keys.
        # TODO - get rid of this. Create it on the fly in the proteins loop.
        readings_by_rep_stage = self._format_protein_readings(protein_readings)

        # TODO - do other optimisations like this
        del(protein_readings)

        raw_readings = self._qc_protein_readings(readings_by_rep_stage)

        for protein, readings in raw_readings.items():
            # logger.info(f"++ PROTEIN: {protein.accession_number}")
            # results[protein] = {
            #     PROTEIN_ABUNDANCES: {
            #         RAW: {},
            #         NORMALISED: {},
            #         IMPUTED: {}
            #     },
            #     PHOSPHORYLATION_ABUNDANCES: {},
            #     METRICS: {}
            # }

            result = {
                PROTEIN_ABUNDANCES: {
                    RAW: {},
                    NORMALISED: {},
                    IMPUTED: {}
                },
                METRICS: {}
            }

            self._calculate_abundances_metrics(
                result,
                project,
                replicates,
                readings,
                column_names,
                sample_stages,
                medians,
                PROTEIN_ABUNDANCES,
                with_bugs
            )

            run_result, _ = RunResult.objects.get_or_create(
                run=run,
                protein=protein
            )

            run_result.protein_result = result

            run_result.save()



    # TODO - put this in the batch process
    def _calculate_batch_q_value_fisher():
        # TODO - this now needs to fetch the anovas from the protein results,
        #   not from the anovas dict.

        relative_log2_readings_by_protein = {}
        anovas = {}

        # fisher_stats = self._calculate_fisher(relative_log2_readings_by_protein, replicates)

        # TODO - should this be outside _proteo?

        # TODO - rename this
        prot_anova_info: dict = {}
        for protein in anovas.keys():
            prot_anova_info[protein] = {}

            prot_anova_info[protein][P_VALUE] = anovas[protein][P_VALUE]

        # TODO - converting to a dataframe seems excessive. Find an alternative.
        # TODO - is this the same as _generate_df?
        prot_anova_info_df = pd.DataFrame(prot_anova_info).T
        prot_anova_info_df[Q_VALUE] = self.p_adjust_bh(prot_anova_info_df[P_VALUE])

        prot_anova_info = prot_anova_info_df.to_dict("index")

        for protein in raw_readings:
            # results[protein][METRICS][LOG2_MEAN][ANOVA] = anovas[protein]

            # ANOVA q values
            q_value = 1

            if protein in prot_anova_info:
                q_value = prot_anova_info[protein][Q_VALUE]

            results[protein][METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = q_value

            # # Fisher
            # fisher = {G_STATISTIC: 1, P_VALUE: 1, FREQUENCY: 1, Q_VALUE: 1}

            # if protein in fisher_stats:
            #     fisher = fisher_stats[protein]

            # results[protein][METRICS][LOG2_MEAN][FISHER_G] = fisher

        return results

    def _calculate_fisher(
        self,
        results,
        replicates,
        phospho=False,
        phospho_ab=False,
        phospho_reg=False,
    ):
        # TODO - remove this comment
        # These calls are always the same in ICR
        # norm_method = LOG2_MEAN
        # raw = False

        time_course_fisher = self._create_results_dataframe(results, phospho, phospho_ab, phospho_reg, replicates)

        time_course_fisher = time_course_fisher.dropna()

        g_stats = []
        p_values = []
        frequencies = []

        for _, row in time_course_fisher.iterrows():
            values = row.values.astype(float)

            # Estimate power spectrum
            freqs, power = periodogram(values)

            if len(power) == 0 or np.all(power == 0):
                g_stat = 0
                p_value = 1.0
                dominant_freq = 0
            else:
                g_stat = np.max(power) / np.sum(power)
                dominant_freq = freqs[np.argmax(power)]
                p_value = 1 - g_stat  # crude approximation

            g_stats.append(g_stat)
            p_values.append(p_value)
            frequencies.append(dominant_freq)

        time_course_fisher[G_STATISTIC] = g_stats
        time_course_fisher[P_VALUE] = p_values
        time_course_fisher[FREQUENCY] = frequencies

        time_course_fisher[Q_VALUE] = self.p_adjust_bh(time_course_fisher[P_VALUE])

        # Return only the Fisher columns
        fisher_cols = [G_STATISTIC, P_VALUE, FREQUENCY, Q_VALUE]
        time_course_fisher = time_course_fisher[fisher_cols]

        return time_course_fisher.to_dict("index")



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



    # TODO - lifted from ICR, rewrite
    # TODO - try to get rid of dataframes if possible
    # TODO - needs a test if to be used, is only used by fisher
    def _create_results_dataframe(self, results, phospho, phospho_ab, phospho_reg, replicates):
        # TODO - raw is always False, ignore all raw code
        # TODO - norm_method is always log2_mean, ignore all norm_method code
        abundance_table = {}
        final_protein = None

        # TODO - maybe this and other similar loops could be changed to use items()?
        #   That way the value variable could be used in each successive loop, e.g.
        #   for protein, p_readings in readings.items():
        #       for replicate_name, rn_readings = p_readings.items():
        #           etc.
        for protein in results:
            prs = results[protein]

            # TODO - just needed for getting replicates I think, could be disposed of
            final_protein = protein

            # TODO - put this 'if' outside the loop? It may be more efficient
            if phospho:

                for mod in prs[PHOSPHORYLATION_ABUNDANCES]:
                    protein_abundances_all = prs[PHOSPHORYLATION_ABUNDANCES][mod][POSITION_ABUNDANCES]

                    mod_key = protein.accession_number + "_" + prs[PHOSPHORYLATION_ABUNDANCES][mod][PHOSPHORYLATION_SITE]

                    if len(protein_abundances_all) != 0:
                        # TODO - tidy this
                        protein_abundances = protein_abundances_all[NORMALISED][LOG2_MEAN]

                        if phospho_ab:
                            if PROTEIN_OSCILLATION_ABUNDANCES in prs[PHOSPHORYLATION_ABUNDANCES][mod]:
                                protein_abundances = prs[PHOSPHORYLATION_ABUNDANCES][mod][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN]
                            else:
                                continue
                        if phospho_reg:
                            if PHOSPHO_REGRESSION in prs[PHOSPHORYLATION_ABUNDANCES][mod]:
                                protein_abundances = prs[PHOSPHORYLATION_ABUNDANCES][mod][PHOSPHO_REGRESSION][LOG2_MEAN]
                            else:
                                continue

                        self._build_phospho_abundance_table(abundance_table, replicates, protein_abundances, mod_key)
            else:
                abundance_table[protein] = {}

                for replicate_name in results[protein]:
                    for stage_name in results[protein][replicate_name].keys():
                        rep_stage_name = f"{replicate_name}_{stage_name}"
                        abundance_table[protein][rep_stage_name] = results[protein][
                            replicate_name
                        ][stage_name]

        time_course_abundance_df = pd.DataFrame(abundance_table)
        time_course_abundance_df = time_course_abundance_df.T

        new_cols = []

        # if phospho:
        #     print("++++++ FOO")
        #     print(results[final_protein][PHOSPHORYLATION_ABUNDANCES].keys())

        #     replicate_names = list(results[final_protein][final_mod].keys())

        #     for stage_name in results[final_protein][final_mod][replicate_names[0]]:
        #         for rn in replicate_names:
        #             new_cols.append(f"{rn}_{stage_name}")
        # else:
        replicate_names = list(results[final_protein].keys())

        for stage_name in results[final_protein][replicate_names[0]]:
            for rn in replicate_names:
                new_cols.append(f"{rn}_{stage_name}")

        # Rearrange the column order so replicates are near eatch other
        try:
            time_course_abundance_df = time_course_abundance_df[new_cols]
        except Exception as e:
            print("++++ DF FAILED")
            print(phospho)
            print(time_course_abundance_df)
            print(e)
            exit()

        return time_course_abundance_df

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
            call_kwargs["replicate_name"] = replicate.name
            results[replicate.name] = func(**call_kwargs)

        return results

    def _tp(self, stage_name, readings):
        """
        Creates a list of each abundance of a stage name across replicates
        """
        res = []

        for replicate_name in readings:
            if stage_name in readings[replicate_name]:
                reading = readings[replicate_name][stage_name]

                if reading is not None:
                    res.append(float(reading))

        return res
    
    # TODO - straight up lifted from ICR, simplify ideally using a library
    def _calculate_ANOVA(self, readings: dict):
        # Defaults if not enough replicate
        # TODO - why these values?
        p_value = 1
        f_statistic = 1

        # TODO - why is it the first replicate?
        first_replicate_name = list(readings.keys())[0]

        # TODO - change to stage_names
        stage_names = []

        try:
            # TODO - change this to use sample_stages?
            for stage_name in readings[first_replicate_name]:
                stage_names.append(self._tp(stage_name, readings))

            # Each entry must have at least two points for f_oneway to work    
            timepoints = [x for x in stage_names if x != [] and len(x) > 1]

            if len(timepoints) > 1:
                # TODO - study f_oneway
                one_way_anova = f_oneway(*timepoints)

                f_statistic = one_way_anova[0].item()
                p_value = one_way_anova[1].item()

                if np.isnan(p_value):
                    p_value = 1

        except Exception as e:
            print("ERROR CALCULATING ANOVA")
            print(e)

        return {
            P_VALUE: p_value,
            F_STATISTICS: f_statistic
        }

    # TODO - tidy this up
    # TODO - lifted from ICR, change names
    # TODO - tested
    def _calculate_metrics(
        self,
        readings: dict,
        readings_averages: dict,
        replicates,
        sample_stages
    ):
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
                if abundance := readings[replicate.name].get(sample_stage.name):
                    # TODO - sometimes the sample_stage.name isn't set. Why?
                    #   Is there something wrong with the populate script?
                    # TODO - why does this add all reps for standard deviation calculation?
                    #   Why isn't it on a per-replicate basis?
                    abundances.append(abundance)

        std = None

        if len(abundances) > 1:
            std = statistics.stdev(abundances)
        # else:
        #     print("NO ABUNDANCES FOR STANDARD DEVIATION")
        #     print(readings)

        # TODO - is this the right thing to do? Go through all of this function
        #   to check the behaviour of its various parts.
        #   Almost all of the warnings output are from this function.
        if not len(abundance_averages_list):
            return metrics

        try:
            residuals_array = np.polyfit(
                range(0, len(abundance_averages)),
                abundance_averages_list,
                2,
                full=True,
            )[1]

            if len(residuals_array) == 0:
                # eg Q9HBL0 {'G2_2': 0.4496, 'G2/M_1': 0.7425, 'M/Early G1': 1.0}
                residuals = 5
            else:
                residuals = np.polyfit(
                range(0, len(abundance_averages)),
                abundance_averages_list,
                2,
                full=True,)[1][0]

            r_squared = self._polyfit(
                range(0, len(abundance_averages)), abundance_averages_list, 2
            )
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
            # TODO - why does it need to be two? Don't we just need the last one?
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
            print("Error in _calculate_metrics")
            print(e)

        return metrics

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

    def _impute(
        # TODO - all these pr types are wrong, and also probably bad variable names
        self,
        readings: dict,
        replicates: QuerySet[Replicate],
        column_names: QuerySet[ColumnName],
    ):
        replicates_by_name: dict = {}
        column_names_by_replicate: dict = {}

        # TODO - is this needed now we no longer use Replicate objects as keys?
        for replicate in replicates:
            replicates_by_name[replicate.name] = replicate
            column_names_by_replicate[replicate.name] = []

        # TODO - and this?
        for column_name in column_names:
            column_names_by_replicate[column_name.replicate.name].append(
                column_name.sample_stage.name
            )

        imputed_readings: dict = {}

        for replicate_name in readings.keys():
            imputed_readings[replicate_name] = {}

            abundances_dict = readings[replicate_name]
            abundances = list(abundances_dict.values())

            stage_names = column_names_by_replicate[replicate_name]

            for idx, stage_name in enumerate(stage_names):
                # Default value, should never be used
                value = 0

                if abundances_dict.get(stage_name) is not None:
                    value = abundances_dict[stage_name]
                else:
                    last = None
                    next = None

                    # TODO - isn't there a better way to iterate?
                    for offset in range(1, len(stage_names)):
                        prev_idx = idx - offset
                        if prev_idx < 0:
                            # Gone before the beginning of the list, give up
                            break

                        prev_stage_name = stage_names[prev_idx]

                        if abundances_dict.get(prev_stage_name) is not None:
                            last = (offset, abundances_dict[prev_stage_name])
                            # last = abundances_dict[prev_stage_name]
                            break

                    for offset in range(1, len(abundances)):
                        # Look forward
                        # TODO - this seems to loop back to the beginning. Is that right?
                        next_idx = (idx + offset) % len(abundances)
                        next_stage_name = stage_names[next_idx]

                        if abundances_dict.get(stage_name) is not None:
                            next = (offset, abundances[next_stage_name])
                            # next = abundances[next_stage_name]
                            break

                    if last and next:
                        # Linear imputation between nearest timepoints
                        # TODO - find out why this calculation
                        # TODO - name variables better
                        d1, a1 = last
                        d2, a2 = next
                        step_height = (a1 - a2) / (d1 + d2)
                        value = d2 * step_height + a2

                # imputed_protein_readings[protein][replicate_name][stage_name] = self._round(float(value))
                # TODO - for some reason ICR rounds to 2, not 4. What to do?
                imputed_readings[replicate_name][stage_name] = round(float(value), 2)

        return imputed_readings

    def _calculate_level_two_normalisation(self, readings: dict, zero_min=False):
        level_two_normalised_readings: dict = {}

        for replicate_name in readings:
            level_two_normalised_readings[replicate_name] = {}

            min_value = 0
            max_value = 0

            abundances = readings[replicate_name]

            abundance_values_non_null = [
                val for val in abundances.values() if val is not None
            ]

            if len(abundance_values_non_null) != 0:
                if not zero_min:
                    min_value = min(abundance_values_non_null)

                max_value = max(abundance_values_non_null)

            for stage_name, abundance in abundances.items():
                denominator = max_value - min_value
                if abundance is None or denominator == 0:
                    abundance_normalised = 0.5
                else:
                    abundance_normalised = (abundance - min_value) / denominator

                level_two_normalised_readings[replicate_name][stage_name] = self._round(
                    abundance_normalised
                )

        return level_two_normalised_readings

    def _calculate_relative_log2_normalisation(self, readings: dict):
        log2_abundances: dict = {}

        for replicate_name in readings:
            log2_abundances[replicate_name] = {}

            for stage_name in readings[replicate_name]:
                log2 = None
                reading = readings[replicate_name][stage_name]

                if reading is not None:
                    log2 = math.log2(reading)

                log2_abundances[replicate_name][stage_name] = log2

        total_abundances = 0
        total_lengths = 0

        log2_normalised_readings: dict = {}

        for replicate_name in log2_abundances:
            for stage_name in log2_abundances[replicate_name]:
                if log2_abundances[replicate_name][stage_name] is not None:
                    total_abundances += log2_abundances[replicate_name][stage_name]
                    total_lengths += 1

            mean = None

            if total_lengths != 0:
                mean = total_abundances / total_lengths
            # TODO - is mean is None then the loop below can be simplified

            log2_normalised_readings = {}

            for replicate_name in readings:
                log2_normalised_readings[replicate_name] = {}

                for stage_name in readings[replicate_name]:
                    normalised_abundance = None

                    if log2_abundances[replicate_name].get(stage_name) is not None and mean is not None:
                        normalised_abundance = self._round(
                            log2_abundances[replicate_name][stage_name] - mean
                        )

                    log2_normalised_readings[replicate_name][
                        stage_name
                    ] = normalised_abundance

        return log2_normalised_readings

    def _calculate_arrest_log2_normalisation(self, readings: dict, project: Project):
        log2_normalised_readings: dict = {}

        # TODO - what should the stage name be?
        # TODO - is ARRESTING_AGENT the wrong name?
        ARRESTING_AGENT = "Nocodozole"

        # TODO - this is a hack, maybe add the field to the Project model?
        if project.name == "ICR":
            ARRESTING_AGENT = "Palbo"

        log2_normalised_readings = {}

        for replicate_name in readings:
            log2_normalised_readings[replicate_name] = {}

            for stage_name in readings[replicate_name]:
                log2_reading = None
                reading = readings[replicate_name][stage_name]

                if readings[replicate_name].get(ARRESTING_AGENT):
                    arrest_reading = readings[replicate_name][ARRESTING_AGENT]

                    if reading is not None and arrest_reading is not None:
                        log2_reading = self._round(math.log2(reading / arrest_reading))

                log2_normalised_readings[replicate_name][stage_name] = log2_reading

        return log2_normalised_readings

    # TODO - why does ICR not need to do QC?
    def _qc_protein_readings(self, all_readings: dict):
        logger.info("Remove any invalid proteins")

        qc_proteins: dict = {}

        for protein in all_readings.keys():
            qc_proteins[protein] = {}

            for replicate_name in all_readings[protein]:
                qc_proteins[protein][replicate_name] = {}

                abundances = all_readings[protein][replicate_name]

                if any(x is not None for x in abundances.values()):
                    qc_proteins[protein][replicate_name] = abundances
                # else:
                #     print(f"+++++ DELETING EMPTY PROTEIN {protein.accession_number}")

        return qc_proteins

    def _format_protein_readings(self, protein_readings: QuerySet[ProteinReading]):
        logger.info(
            "Converting protein_readings QuerySet into dict by protein, replicate and stage name"
        )

        readings_by_rep_stage: dict = {}

        protein_no = 0

        for protein_reading in protein_readings:
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Formatting protein for {protein_no}, {protein_reading.protein.accession_number}",
            )

            # TODO - what about Nones? Will there be any here? Check the import script.
            reading = protein_reading.reading

            protein = protein_reading.protein
            replicate_name = protein_reading.column_name.replicate.name
            stage_name = protein_reading.column_name.sample_stage.name

            if not readings_by_rep_stage.get(protein):
                readings_by_rep_stage[protein] = {}

            if not readings_by_rep_stage[protein].get(replicate_name):
                readings_by_rep_stage[protein][replicate_name] = {}

            readings_by_rep_stage[protein][replicate_name][stage_name] = reading

        return readings_by_rep_stage

    def _format_phospho_readings(self, phospho_readings: QuerySet[PhosphoReading]):
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

    def _get_phosphosites(self, phosphos: QuerySet[Phospho]):
        logger.info(
            "Get all the phosphites"
        )
        phosphosites: dict = {}

        mod_no = 0

        for phospho in phosphos:
            mod_no += 1
            self._count_logger(
                mod_no,
                10000,
                f"Formatting phosphites for {phospho.mod}",
            )

            phosphosites[phospho.mod] = phospho.phosphosite

        return phosphosites

    def _calculate_first_level_normalisation(self, readings: dict, medians):
        normalised_readings: dict = {}

        for replicate_name in readings:
            normalised_readings[replicate_name] = {}

            for stage_name in readings[replicate_name]:
                normalised_reading = None

                reading = readings[replicate_name][stage_name]

                if reading is not None:
                    median = medians[replicate_name][stage_name]

                    # TODO - need to round this?
                    normalised_reading = reading / median

                normalised_readings[replicate_name][stage_name] = normalised_reading

        return normalised_readings

    # calculateAverageRepAbundance
    def _calculate_means(
        self,
        reading: dict,
        with_bugs: bool,
        imputed: bool = False,
    ):
        means: dict = {}

        abundances: dict = {}

        for replicate_name in reading:
            for stage_name in reading[replicate_name]:
                if not abundances.get(stage_name):
                    abundances[stage_name] = []

                if with_bugs and not imputed:
                    # We throw away the second reading
                    # TODO - how will this behave towards None?
                    if len(abundances[stage_name]) == 1:
                        continue

                if reading[replicate_name][stage_name] is not None:
                    abundances[stage_name].append(reading[replicate_name][stage_name])

            means = {}

            for stage_name in abundances:
                abundance = abundances[stage_name]

                if len(abundance):
                    mean = sum(abundance) / len(abundance)
                    means[stage_name] = self._round(mean)
                else:
                    # TODO - is this the right thing to do?
                    means[stage_name] = None

        return means

    def _calculate_phospho_medians(
        self,
        phospho_readings: dict,
        run
    ):
        # TODO - get rid of _format_phospho_readings
        #   Correct output:
        #   {'One': {'Palbo': 83.6, 'Late G1_1': 90.7, 'G1/S': 73.4, 'S': 76.6, 'S/G2': 63.8, 'G2_2': 91.75, 'G2/M_1': 107.0, 'M/Early G1': 181.8}, 'Two': {'Palbo': 103.95, 'Late G1_1': 71.2, 'G1/S': 122.8, 'S': 77.9, 'S/G2': 89.3, 'G2_2': 116.0, 'G2/M_1': 122.9, 'M/Early G1': 122.7}}

        readings = self._format_phospho_readings(phospho_readings)

        medians = {}

        for protein in readings.keys():
            for mod in readings[protein].keys():
                for replicate_name in readings[protein][mod].keys():
                    if not medians.get(replicate_name):
                        medians[replicate_name] = {}

                    for column_name in readings[protein][mod][replicate_name].keys():
                        if not medians[replicate_name].get(column_name):
                            medians[replicate_name][column_name] = []

                        reading = readings[protein][mod][replicate_name][column_name]

                        if reading is None:
                            continue

                        medians[replicate_name][column_name].append(reading)

        for replicate_name in medians.keys():
            for column_name in medians[replicate_name].keys():
                median = statistics.median(medians[replicate_name][column_name])

                medians[replicate_name][column_name] = median

        run, _ = Run.objects.get_or_create(
            project=run.project, with_bugs=run.with_bugs, limit_proteins=False
        )

        run.phospho_medians = json.dumps(medians)

        run.save()

        return medians


    def _calculate_replicate_stage_name_medians(
        self,
        replicate_name: str,
        protein_readings: QuerySet[ProteinReading],
        column_names: QuerySet[ColumnName],
    ):
        stage_name_medians = {}

        column_names_by_replicate = column_names.filter(replicate__name=replicate_name)

        for column_name in column_names_by_replicate:
            readings = []

            protein_readings_by_column = protein_readings.filter(
                column_name=column_name
            )

            for protein_reading in protein_readings_by_column:
                # TODO - what to do about None values?
                if protein_reading.reading is not None:
                    readings.append(protein_reading.reading)

            if len(readings) == 0:
                raise Exception(
                    "Can't create median with no abundances for protein {protein_reading.protein.accession_number}"
                )

            median = statistics.median(readings)

            stage_name_medians[column_name.sample_stage.name] = median

        return stage_name_medians

    # TODO - think this could be tidied
    # TODO - tested
    def _generate_xs_ys(self, readings, replicates):
        # TODO - why do we use the last replicate? It's in ICR
        final_replicate = list(replicates)[-1]

        stage_names_map = {}

        # TODO - change this to use samples_stages list
        for i, stage_name in enumerate(readings[final_replicate.name].keys()):
            stage_names_map[stage_name] = i

        x = []
        for stage_name in readings.get(final_replicate.name, {}):
            x.append(stage_names_map[stage_name])
        x.sort()

        y = []
        for stage_name in stage_names_map:
            value = readings.get(final_replicate.name, {}).get(stage_name)
            if value is not None:
                y.append(value)

        return x, y, stage_names_map
    
    # TODO - lifted from ICR, change
    def _getConsensusKinasePred(self, uniprot_accession, phosphosite):
        """
        Predicts the kinase most likely to phosphorylate a phosphorylation site 
        based on the consensus approach.
        """
        phospho_kinases_class = {}
        peptide_seq = self._getPeptideAligned(uniprot_accession, phosphosite)

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
            phospho_kinases_class = {"accession":uniprot_accession,"site":phosphosite, "peptide_seq": peptide_seq, "kinase_motif_match": motif_matches}
        if matches == 0:
            total_not_matched+=1
            phospho_kinases_class = {"accession":uniprot_accession, "site":phosphosite, "peptide_seq": peptide_seq, "kinase_motif_match": ["-"]}

        return phospho_kinases_class

    # TODO - lifted from ICR
    def _getPeptideAligned(self, uniprot_accession, phosphosite):
        """
        Checks if the phosphosite is centered in the peptide sequence and it aligns it if not.
        """
        aa = phosphosite[0]
        position = int(phosphosite[1::]) -1
        site =  str(position+1)

        if 'phospho' in index_protein_names[uniprot_accession]:
            if site in index_protein_names[uniprot_accession]['phospho']:
                if 'peptide_seq' in index_protein_names[uniprot_accession]['phospho'][site]:
                    peptide_seq = index_protein_names[uniprot_accession]['phospho'][site]['peptide_seq']
                else:
                    peptide_seq = self._getPeptideSequence(uniprot_accession, phosphosite)

        phospho_alignment = ""

        #Discard invalid inputs
        if len(peptide_seq) == 0:           
            peptide_seq = index_protein_names[uniprot_accession]['phospho'][site]['Peptide']
            phospho_alignment = peptide_seq[5:16]
        # Middle of the protein sequence
        elif len(peptide_seq) == 11:
            if aa == peptide_seq[5]:
                phospho_alignment = peptide_seq
            else:
                # Site not in the middle of seq
                peptide_seq_new = index_protein_names[uniprot_accession]['phospho'][site]['Peptide']
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
    def _getPeptideSequence(self, uniprot_accession, phosphosite):
        """
        Creates a peptide sequence which is a substring of the original protein sequence.
        +5/-5
        amino acids from the phosphorylation site.
        """
        sequence = self._getProteinSeq(uniprot_accession)

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
    def _getProteinSeq(self, uniprot_accession):
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
    def _add_protein_oscillations(self, results, replicates, sample_stages, with_bugs):
        # TODO - check all logging statements for similarity to ICR
        logger.info("Adding Protein Oscillation Normalised Abundances")

        self._generate_protein_metrics(results, replicates, sample_stages, with_bugs)

        # Fisher G Statistic - Phospho
        # TODO - figure out how to get this working for protein oscillations
        # time_course_fisher_dict = self._calculate_fisher(results, phospho = True, phospho_ab = True)

        # return

        # Corrected q values - Phospho
        # 1) Create a dataframe with the desired Protein-Phospho info
        prot_phospho_info = {}
        for protein in results:
            if len(results[protein][PHOSPHORYLATION_ABUNDANCES]) != 0:
                for site in results[protein][PHOSPHORYLATION_ABUNDANCES]:
                    if 'protein_oscillation_abundances' in results[protein][PHOSPHORYLATION_ABUNDANCES][site]:
                        phospho_key = protein.accession_number + "_" + results[protein][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]
                        if phospho_key not in prot_phospho_info:
                            prot_phospho_info[phospho_key] = {}
                        prot_phospho_info[phospho_key]['p_value'] = results[protein][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][P_VALUE]

        if not prot_phospho_info:
            # TODO - make this an exception? Shouldn't happen for large data sets
            print("++++ NO PROT_PHOSPHO_INFO")
            exit()
            return results

        prot_phospho_info = self._generate_df(prot_phospho_info)

        # 4) Add Protein-Phospho info in combined_time_course_info dictionary
        # TODO - tidy up, two loops not necessary?
        for protein in results:
            if len(results[protein][PHOSPHORYLATION_ABUNDANCES]) != 0 and len(results[protein][PROTEIN_ABUNDANCES][RAW]) != 0:
                for site in results[protein][PHOSPHORYLATION_ABUNDANCES]:
                    if PROTEIN_OSCILLATION_ABUNDANCES in results[protein][PHOSPHORYLATION_ABUNDANCES][site]:
                        site_key = protein.accession_number + "_" + results[protein][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]

                        # ANOVA q values
                        q_value = 1

                        if site_key in prot_phospho_info:
                            q_value = prot_phospho_info[site_key]['q_value']
                        
                        results[protein][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = q_value

                        # # Fisher
                        # TODO - tidy this and others like it
                        # if site_key in time_course_fisher_dict:
                        #     results[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS]["Fisher_G"] = time_course_fisher_dict[site_key]
                        # else:
                        #     results[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS]["Fisher_G"] = {'G_statistic': 1, 'p_value': 1, 'frequency': 1, 'q_value': 1}

        return results





    # TODO - lifted from ICR
    # addPhosphoRegression
    def _add_phospho_regression(self, results, replicates, sample_stages, with_bugs):
        # TODO - check and revise all comments
        # Add the Regressed Phospho Normalised Abundances
        """
        Normalise the phospho abundance on the protein abundance
        Calculates and adds all the Regressed Phospho Abundance and their metrics for each phosphosite.
        Phospho = Dependent = Y
        Protein = Independent = X
        Linear Model => y = ax + b
        Residuals = Y - Y_predict 
        """
        logger.info("Adding Phospho Normalised on Protein Abundances - Regression")

        self._generate_phospho_regression_metrics(results, replicates, sample_stages, with_bugs)

        # # Fisher G Statistic - Phospho
        # time_course_fisher_dict = self._calculate_fisher(combined_time_course_info, phospho = True, phospho_ab = False, phospho_reg = True)

        # Corrected q values - Phospho
        # 1) Create a dataframe with the desired regression info
        regression_info = {}

        for protein in results:
            if len(results[protein][PHOSPHORYLATION_ABUNDANCES]) != 0:
                for site in results[protein][PHOSPHORYLATION_ABUNDANCES]:
                    if 'phospho_regression' in results[protein][PHOSPHORYLATION_ABUNDANCES][site]:
                        # TODO - a hack
                        phospho_key = protein.accession_number + "_" + results[protein][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]
                        if phospho_key not in regression_info:
                            regression_info[phospho_key] = {}

                        if results[protein][PHOSPHORYLATION_ABUNDANCES][site]['phospho_regression']['log2_mean'].get(METRICS):
                            # TODO - why would it not have metrics? Lack of data maybe?
                            regression_info[phospho_key]['p_value'] = results[protein][PHOSPHORYLATION_ABUNDANCES][site]['phospho_regression']['log2_mean'][METRICS]['ANOVA']['p_value']

        regression_info = self._generate_df(regression_info)

        # 4) Add Regression info in combined_time_course_info dictionary
        for protein in results:
            if len(results[protein][PHOSPHORYLATION_ABUNDANCES]) != 0 and len(results[protein][PROTEIN_ABUNDANCES][RAW]) != 0:
                for phosphosite in results[protein][PHOSPHORYLATION_ABUNDANCES]:
                    if PHOSPHO_REGRESSION in results[protein][PHOSPHORYLATION_ABUNDANCES][phosphosite]:
                        site_key = protein.accession_number + "_" + results[protein][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHORYLATION_SITE]

                        # ANOVA q values
                        q_value = 1

                        if site_key in regression_info:
                            q_value = regression_info[site_key]['q_value']

                        if results[protein][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN].get(METRICS):
                            # TODO - why no metrics? Not generated due to lack of data?
                            results[protein][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = q_value

                        # # Fisher
                        # phospho_regression[LOG2_MEAN][METRICS]["Fisher_G"] = {}
                        # if site_key in time_course_fisher_dict:
                        #     combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS]["Fisher_G"] = time_course_fisher_dict[site_key]
                        # else:
                        #     combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
                        #                 phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS]["Fisher_G"] = {'G_statistic': 1, 'p_value': 1, 'frequency': 1, 'q_value': 1}

        return results

    # calculateProteinOscillationAbundances
    # TODO - tested
    def _calculate_protein_oscillation(self, results_single, norm_method, phosphosite, replicates):
        phospho_oscillations = {}

        rspannm = results_single[PROTEIN_ABUNDANCES][NORMALISED][norm_method]
        rspappann = results_single[PHOSPHORYLATION_ABUNDANCES][
                        phosphosite][POSITION_ABUNDANCES][NORMALISED][norm_method]

        for replicate in replicates:
            if replicate.name in rspannm:
                if (replicate.name in rspappann):
                    protein_oscillation_abundances = {}
                    protein_normed = rspannm[replicate.name]

                    phospho_normed = rspappann[replicate.name]

                    try:
                        for stage_name in protein_normed:
                            if not protein_normed[stage_name] or not phospho_normed.get(stage_name):
                                continue

                            protein_oscillation_abundances[stage_name] = (
                                phospho_normed[stage_name]
                                - protein_normed[stage_name]
                            )
                    except Exception as e:
                        print("+++++ _calculate_protein_oscillation error")
                        print(protein_normed)
                        print(phospho_normed)
                        print(e)
                        exit()

                    phospho_oscillations[replicate.name] = protein_oscillation_abundances
    
        return phospho_oscillations

    def _calculate_abundances_metrics(
        self,
        result,
        project,
        replicates,
        readings,
        column_names,
        sample_stages,
        medians,
        location,
        with_bugs
    ):
        # firstLevelNormalisationProteomics
        # firstLevelNormalisationPhospho
        normalised_readings = self._calculate_first_level_normalisation(
            readings, medians
        )

        # calclog2PalboNormalisation
        arrest_readings = self._calculate_arrest_log2_normalisation(
            normalised_readings, project
        )

        # calclog2RelativeAbundance
        log2_readings = (
            self._calculate_relative_log2_normalisation(normalised_readings)
        )

        # normaliseData
        min_max_readings = self._calculate_level_two_normalisation(
            log2_readings
        )

        # normaliseData
        zero_max_readings = self._calculate_level_two_normalisation(
            normalised_readings, True
        )


        imputed_readings = self._impute(
            min_max_readings, replicates, column_names
        )

        raw_averages = (
            self._calculate_means(
                readings, with_bugs
            )
        )

        normalised_averages = (
            self._calculate_means(
                normalised_readings, with_bugs
            )
        )

        min_max_averages = (
            self._calculate_means(
                min_max_readings, with_bugs
            )
        )

        zero_max_averages = (
            self._calculate_means(
                zero_max_readings, with_bugs
            )
        )

        log2_averages = (
            self._calculate_means(
                log2_readings, with_bugs
            )
        )

        arrest_averages = (
            self._calculate_means(
                arrest_readings, with_bugs
            )
        )

        imputed_averages = (
            self._calculate_means(
                imputed_readings, imputed=True, with_bugs = with_bugs
            )
        )

        log2_mean_metrics = self._calculate_metrics(
            log2_readings,
            log2_averages,
            replicates,
            sample_stages
        )

        zero_max_mean_metrics = self._calculate_metrics(
            zero_max_readings,
            zero_max_averages,
            replicates,
            sample_stages
        )

        # # TODO - move this out, it is for multiple proteins
        anova = self._calculate_ANOVA(log2_readings)

        # relative_log2_readings_by_protein[
        #     protein
        # ] = log2_readings

        result[location][RAW] = readings
        result[location][RAW][ABUNDANCE_AVERAGE] = raw_averages
        # TODO - confirm the output later calculations are as they should be after this
        result[location][NORMALISED][LOG2_MEAN] = copy.deepcopy(log2_readings)
        result[location][NORMALISED][LOG2_MEAN][ABUNDANCE_AVERAGE] = log2_averages
        result[location][NORMALISED][MIN_MAX] = min_max_readings
        result[location][NORMALISED][MIN_MAX][ABUNDANCE_AVERAGE] = min_max_averages
        result[location][NORMALISED][MEDIAN] = normalised_readings
        result[location][NORMALISED][MEDIAN][ABUNDANCE_AVERAGE] = normalised_averages
        result[location][NORMALISED][ZERO_MAX] = zero_max_readings
        result[location][NORMALISED][ZERO_MAX][ABUNDANCE_AVERAGE] = zero_max_averages
        result[location][NORMALISED][LOG2_ARREST] = arrest_readings
        result[location][NORMALISED][LOG2_ARREST][ABUNDANCE_AVERAGE] = arrest_averages
        result[location][IMPUTED] = imputed_readings
        result[location][IMPUTED][ABUNDANCE_AVERAGE] = imputed_averages
        result[METRICS][LOG2_MEAN] = log2_mean_metrics
        result[METRICS][LOG2_MEAN][ANOVA] = anova
        result[METRICS][ZERO_MAX] = zero_max_mean_metrics



    # TODO - tested
    def _generate_phospho_regression_metrics(self, results, replicates, sample_stages, with_bugs):
        stages = [s.name for s in sample_stages]

        for protein in results:
            phosphoryl_abundances = results[protein][PHOSPHORYLATION_ABUNDANCES]
            protein_abundances = results[protein][PROTEIN_ABUNDANCES]

            try:
                if len(phosphoryl_abundances) == 0 or len(protein_abundances[RAW]) == 0:
                    continue

                for phosphosite in phosphoryl_abundances:
                    # TODO - not in the original
                    pappan = phosphoryl_abundances[
                        phosphosite][POSITION_ABUNDANCES][NORMALISED]

                    rep_proteos = []
                    rep_phosphos = []

                    for replicate in replicates:
                        if not pappan[LOG2_MEAN].get(replicate.name):
                            continue

                        rep_proteo = protein_abundances[NORMALISED][LOG2_MEAN][replicate.name]
                        rep_phospho = pappan[LOG2_MEAN][replicate.name]

                        if len(rep_proteo) != len(sample_stages) or len(rep_phospho) != len(sample_stages):
                            continue

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
                    phospho_Y_list = [val for r_phos in rep_phosphos for val in r_phos.values()]
                    protein_X_list = [val for r_prot in rep_proteos for val in r_prot.values()]

                    # converting list to array
                    prot_x = np.asarray(protein_X_list).reshape(len(protein_X_list), 1)
                    phospho_y = np.asarray(phospho_Y_list).reshape(len(phospho_Y_list), 1)

                    try:
                        # Fit the linear model
                        if None in prot_x or None in phospho_y:
                            continue

                        model = linear_model.LinearRegression().fit(prot_x, phospho_y)
                    except Exception as e:
                        print("++++++ LINEAR REGRESSION")
                        print(prot_x)
                        print(phospho_y)
                        print(e)

                    # Predict new Phospho Values
                    y_pred = model.predict(prot_x)

                    # Calculate Residuals
                    residuals = (phospho_y - y_pred)

                    # Create new regressed phospho abundances dictionaries
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

                    # Add metrics

                    # Calculate the protein - phospho vector correlation
                    phospho_regression[ZERO_MAX][METRICS] = {}

                    # [0] to get the correlation coefficient, [1] = p-value
                    phospho_regression[ZERO_MAX][METRICS][PROTEIN_PHOSPHO_CORRELATION] = stats.pearsonr(protein_X_list, phospho_Y_list)[0]

                    # curve fold change phosphorylation/curve fold change protein for 0-max
                    phospho_regression[ZERO_MAX][METRICS][PHOSPHO_PROTEIN_CFC_RATIO] = results[protein][PHOSPHORYLATION_ABUNDANCES][
                        phosphosite][METRICS][ZERO_MAX][CURVE_FOLD_CHANGE] / results[protein][METRICS][ZERO_MAX][CURVE_FOLD_CHANGE]

                    # ANOVA
                    phospho_regression[LOG2_MEAN][METRICS] = {}
                    phospho_regression[LOG2_MEAN][METRICS][ANOVA] = self._calculate_ANOVA(phospho_regression[LOG2_MEAN])
            except Exception as e:
                print("+++++ ERROR")
                print(e)
                exit()

        return results

    # TODO - tested
    def _generate_protein_metrics(self, results, replicates, sample_stages, with_bugs):
        for protein in results:
            # If we have info both in protein and in phospho level
            prpa = results[protein][PHOSPHORYLATION_ABUNDANCES]

            if len(prpa) != 0 and len(results[protein][PROTEIN_ABUNDANCES][RAW]) != 0:
                for mod in prpa:
                    prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES] = {ZERO_MAX:{}, LOG2_MEAN:{}}
                    prpampoa = prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES]

                    for norm_method in [ZERO_MAX, LOG2_MEAN]:
                        protein_oscillation_abundances = prpampoa[norm_method]

                        # TODO - should this be passing phosphosite, not mod?
                        phospho_oscillations = self._calculate_protein_oscillation(
                            results[protein], norm_method, mod, replicates
                        )

                        for rep in replicates:
                            if rep.name in phospho_oscillations:
                                protein_oscillation_abundances[rep.name] = phospho_oscillations[rep.name]

                        # TODO - split this out into its own variable? That's done with other average
                        protein_oscillation_abundances[ABUNDANCE_AVERAGE] = self._calculate_means(
                            protein_oscillation_abundances, imputed = False, with_bugs = with_bugs
                        )

                    # if we have info in Protein Oscillation Normalised Abundances
                    # TODO - check all comments to make sure they're non-ICR
                    if len(prpampoa[LOG2_MEAN]) > 1:
                        for norm_method in [ZERO_MAX, LOG2_MEAN]:
                            protein_oscillation_abundances = prpampoa[norm_method]

                            # Metrics
                            prpampoa[norm_method][METRICS] = self._calculate_metrics(
                                protein_oscillation_abundances,
                                protein_oscillation_abundances[ABUNDANCE_AVERAGE],
                                replicates,
                                sample_stages
                            )

                        anovas = self._calculate_ANOVA(prpampoa[LOG2_MEAN])

                        prpampoa[LOG2_MEAN][METRICS][ANOVA] = anovas

        return results

    def _merge_phospho_with_proteo(self, results, phospho_results):
        # Not all phospho proteins are in the protein results, add them if required
        for protein in phospho_results:
            # TODO - this can be tidied, has duplication
            # TODO - lots of code should be put in functions and tested
            if protein in results:
                results[protein][
                    PHOSPHORYLATION_ABUNDANCES
                ] = phospho_results[protein]
            else:
                # TODO - put these in as needed
                # if protein in index_protein_names:
                #     gene_name = index_protein_names[protein.accession_number]["gene_name"]
                #     protein_name = index_protein_names[protein.accession_number][
                #         "protein_name"
                #     ]
                # else:
                #     gene_name, protein_name = getProteinInfo(protein.accession_number)

                results[protein] = {
                    PROTEIN_ABUNDANCES: {RAW: {}, NORMALISED: {}, IMPUTED: {}},
                    PHOSPHORYLATION_ABUNDANCES: phospho_results[protein],
                }
                # "gene_name": gene_name,
                # "protein_name": protein_name,

        return results

    def _build_phospho_abundance_table(abundance_table, replicates, protein_abundances, mod_key):
        # TODO - tidy this
        print("+++++ BUILD")
        print(abundance_table)
        print(replicates)
        print(protein_abundances)
        print(mod_key)


        for rep in replicates:
            # TODO - is the "abundance_" bit really needed?
            abundance_rep = f"abundance_{rep.name}"
            if mod_key not in abundance_table:
                abundance_table[mod_key] = {}

            # TODO - not in original, why here?
            if not protein_abundances.get(abundance_rep):
                continue

            # TODO - change timepoint to stage_name
            # TODO - use sample_stages
            for timepoint in protein_abundances[abundance_rep]:
                rep = "_".join(abundance_rep.split("_", 1)[1].split("_")[:2])
                rep_timepoint = rep + "_" + timepoint
                abundance_table[mod_key][rep_timepoint] = protein_abundances[abundance_rep][timepoint]

        print("+++++ BUILD")
        print(abundance_table)
        exit()

        return abundance_table

    # TODO - tested
    def _generate_df(self, info):
        info_df = pd.DataFrame(info)
        info_df = info_df.T

        # 2) Regression ANOVA q values
        info_df['q_value'] = self.p_adjust_bh(info_df['p_value'])

        # 3) Turn dataframe into a dictionary
        return info_df.to_dict('index')

    def _generate_kinase(self, raw_readings, results):
        # Kinase Consensus Prediction
        # TODO - lifted, change
        for protein in raw_readings:
            for mod in raw_readings[protein]:
                # only for certain phosphorylation sites
                # TODO - why? What does the lack of dash mean?
                if mod.find("-") == -1:
                    # Add PepTools Phospho annotations
                    mod_result = results[protein][mod]

                    if 'phospho' in index_protein_names[protein.accession_number] and mod in index_protein_names[protein.accession_number]['phospho']:
                        mod_result['PepTools_annotations'] = index_protein_names[protein.accession_number]['phospho'][mod]

                    mod_result['kinase_prediction'] = {}

                    phospho_kinases_class = self._getConsensusKinasePred(protein.accession_number, mod_result)
                    mod_result['kinase_prediction']['peptide_seq'] = phospho_kinases_class['peptide_seq']
                    mod_result['kinase_prediction']['consenus_motif_match'] = phospho_kinases_class['kinase_motif_match']


    def _calculate_protein_medians(self, run, replicates, protein_readings, column_names, with_bugs):
        # TODO - rename 'medians' to something more informative?
        # TODO - does this need to be by replicate? Why not just all columns at once?
        # TODO - is _all_replicates really useful?
        medians = self._all_replicates(
            func=self._calculate_replicate_stage_name_medians,
            replicates=replicates,
            protein_readings=protein_readings,
            column_names=column_names,
        )

        if with_bugs:
            medians["One"] = medians["Two"]

        run, _ = Run.objects.get_or_create(
            project=run.project, with_bugs=run.with_bugs, limit_proteins=False
        )

        run.protein_medians = json.dumps(medians)

        return medians

    def _fetch_protein_medians(self, run):
        # Load the proteo medians for this project.
        # N.B. THIS LOADS THE MEDIANS FOR ALL PROTEINS!
        #   Not for limit_proteins, that has no use.
        unlimited_run = Run.objects.get(
            project=run.project, with_bugs=run.with_bugs, limit_proteins=False
        )

        medians = unlimited_run.protein_medians

        if not medians:
            raise Exception(f"No medians created yet for {unlimited_run}")

        return json.loads(medians)

    def _fetch_phospho_medians(self, run):
        # Load the phospho medians for this project.
        # N.B. THIS LOADS THE MEDIANS FOR ALL PHOSPHOS!
        #   Not for limit_proteins, that has no use.
        unlimited_run = Run.objects.get(
            project=run.project, with_bugs=run.with_bugs, limit_proteins=False
        )

        medians = unlimited_run.phospho_medians

        if not medians:
            raise Exception(f"No medians created yet for {unlimited_run}")

        return json.loads(medians)


    def _save_data(self, results, file, results_based = True):
        def convert_np(obj):
            if isinstance(obj, np.generic):
                return obj.item()
            if isinstance(obj, float) and math.isnan(obj):
                return None

            return str(obj)
        
        stripped = {}

        if results_based:
            for protein in results:
                stripped[protein.accession_number] = results[protein]
        else:
            stripped = results

        with open(f"output/{file}", "w") as outfile:
            json.dump(stripped, outfile, default=convert_np)

    def _dump(self, obj):
        print(json.dumps(obj, default=lambda o: o.item() if isinstance(o, np.generic) else str(o)))

    def _count_logger(self, i: int, step: int, output: str):
        if i % step == 0:
            logger.info(output)

    def _round(self, value):
        return round(value, 4)
