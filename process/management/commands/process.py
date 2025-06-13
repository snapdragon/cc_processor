import json
import logging
import math
import statistics
import copy

import numpy as np
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from django.db.models.query import QuerySet
from scipy.signal import periodogram
from scipy.stats import f_oneway, moment
from sklearn.metrics import r2_score

from process.models import ColumnName, Project, Protein, ProteinReading, Replicate, PhosphoReading

# TODO - make this configurable by flag?
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

# TODO - move constants elsewhere
FOCUS_PROTEIN_ACCESSION_NUMBER = "Q09666"
RAW = "raw"
METRICS = "metrics"
LOG2_MEAN = "log2_mean"
ZERO_MAX = "0-max"
ANOVA = "ANOVA"
Q_VALUE = "q value"
FISHER_G = "Fisher G"
PROTEIN_ABUNDANCES = "protein_abundances"
ABUNDANCE_AVERAGE = "average"
NORMALISED = "normalised"
MEDIAN = "median"
MIN_MAX = "min-max"
# TODO - change this name, it's specific to ICR
LOG2_PALBO = "log2 palbo"
IMPUTED = "imputed"

class Command(BaseCommand):
    help = "Processes all proteins for a given project"

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

    def handle(self, *args, **options):
        project_name = options["project"]
        with_bugs = options["with_bugs"]

        if with_bugs and project_name != "ICR":
            raise CommandError("Only an ICR project can run --with-bugs")

        logger.info(f"Processing for project {project_name}, with bugs {with_bugs}")

        project = Project.objects.get(name=project_name)
        replicates = Replicate.objects.filter(project=project)
        column_names = ColumnName.objects.filter(replicate__project=project)
        protein_readings = ProteinReading.objects.filter(
            column_name__replicate__project=project
        )

        # self._proteo(project, replicates, protein_readings, column_names, with_bugs)

        phospho_readings = PhosphoReading.objects.filter(phospho__protein__project=project)
        # phospho_readings = PhosphoReading.objects.filter(phospho__protein__project=project)[:100]

        self._phospho(project, replicates, phospho_readings, column_names, with_bugs)



    def _phospho(self, project, replicates, phospho_readings, column_names, with_bugs: bool):
        logger.info("Processing phosphoproteome")

        # TODO - convert this in to a DB query to save having to run it manually?
        readings_by_rep_stage = self._format_phospho_readings(phospho_readings)

        # raw_readings = self._qc_phospho_readings(readings_by_rep_stage)

        raw_readings = readings_by_rep_stage

        medians = self._calculate_phospho_medians(raw_readings)

        # print("++++ MEDIANS")
        # print(medians)
        # exit()

        num_proteins = 0

        results = {}

        for protein in raw_readings.keys():
            for mod, readings in raw_readings[protein].items():
                num_proteins += 1

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
                min_max_normalised_readings = self._calculate_level_two_normalisation(
                    normalised_readings
                )

                # normaliseData
                zero_max_normalised_readings = self._calculate_level_two_normalisation(
                    normalised_readings, True
                )

                imputed_readings = self._impute(
                    min_max_normalised_readings, replicates, column_names
                )

                raw_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        readings, with_bugs
                    )
                )

                min_max_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        min_max_normalised_readings, with_bugs
                    )
                )

                zero_max_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        zero_max_normalised_readings, with_bugs
                    )
                )

                log2_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        log2_readings, with_bugs
                    )
                )

                # TODO - why is this called median when it's from normalised?
                median_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        normalised_readings, with_bugs
                    )
                )

                arrest_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        arrest_readings, with_bugs
                    )
                )

                imputed_averages = (
                    self._calculate_means_across_replicates_by_stage(
                        imputed_readings, imputed=True, with_bugs = with_bugs
                    )
                )

                # print("+++++ IMPUTED AVERAGES")
                # print(protein)
                # print(mod)
                # print(imputed_averages)
                # exit()

        # Metrics
        # time_course_phospho_full = self._addPhosphoMetrics(time_course_phospho_full)



    def _proteo(
        self, project, replicates, protein_readings, column_names, with_bugs: bool
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

        # TODO - make this a flag
        FOCUS_PROTEIN = Protein.objects.get(
            accession_number=FOCUS_PROTEIN_ACCESSION_NUMBER, project__name=project.name
        )

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
            r2_medians = {}

            # TODO - is this still needed now we use replicate names instead of replicates?
            for stage_name in medians["Two"].keys():
                r2_medians[stage_name] = medians["Two"][stage_name]

            for stage_name in medians["One"].keys():
                medians["One"][stage_name] = r2_medians[stage_name]

        # TODO - is this used for anything?
        # means_across_replicates_by_stage = (
        #     self._calculate_means_across_replicates_by_stage(
        #         protein_readings, with_bugs
        #     )
        # )

        # N.B. protein_readings_by_rep_stage is not the same structure as protein_readings.
        #   protein_readings is just a list of ProteinReading objects. normalised_protein_readings
        #   is a dict with Protein object keys. The values is a dict of replicate name keys
        #   with a dict of sample stage names and abundances as keys.
        # TODO - convert this in to a DB query to save having to run it manually?
        readings_by_rep_stage = self._format_protein_readings(protein_readings)

        raw_readings = self._qc_protein_readings(readings_by_rep_stage)

        num_proteins = 0

        relative_log2_readings_by_protein = {}
        anovas = {}
        results = {}

        for protein, readings in raw_readings.items():
            # logger.info(f"++ PROTEIN: {protein.accession_number}")
            num_proteins += 1

            results[protein] = {
                PROTEIN_ABUNDANCES: {
                    RAW: {},
                    NORMALISED: {},
                    IMPUTED: {}
                },
                METRICS: {}
            }

            # TODO - check whether all means calculations need with-bugs
            # TODO - change all these to 'raw_averages' and the like
            raw_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    readings, with_bugs
                )
            )

            # firstLevelNormalisationProteomics
            normalised_readings = self._calculate_first_level_normalisation(
                readings, medians
            )

            normalised_means_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    normalised_readings, with_bugs
                )
            )

            # calclog2PalboNormalisation
            arrest_log2_normalised_readings = self._calculate_arrest_log2_normalisation(
                normalised_readings, project
            )

            arrest_log2_normalised_means_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    arrest_log2_normalised_readings, with_bugs
                )
            )

            # calclog2RelativeAbundance
            relative_log2_normalised_readings = (
                self._calculate_relative_log2_normalisation(normalised_readings)
            )

            relative_log2_normalised_means_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    relative_log2_normalised_readings, with_bugs
                )
            )

            # normaliseData
            # TODO - rename these to min-max and zero-max. Also they're back to front.
            level_two_normalised_readings = self._calculate_level_two_normalisation(
                relative_log2_normalised_readings
            )

            level_two_normalised_means_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    level_two_normalised_readings, with_bugs
                )
            )

            min_max_normalised_readings = self._calculate_level_two_normalisation(
                normalised_readings, True
            )

            min_max_normalised_means_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    min_max_normalised_readings, with_bugs
                )
            )

            imputed_readings = self._impute(
                level_two_normalised_readings, replicates, column_names
            )

            imputed_across_replicates_by_stage = (
                self._calculate_means_across_replicates_by_stage(
                    imputed_readings, with_bugs, imputed=True
                )
            )

            log2_mean_metrics = self._calculate_metrics(
                relative_log2_normalised_readings,
                relative_log2_normalised_means_across_replicates_by_stage,
            )

            zero_max_mean_metrics = self._calculate_metrics(
                min_max_normalised_readings,
                min_max_normalised_means_across_replicates_by_stage,
            )

            anovas[protein] = self._calcANOVA(relative_log2_normalised_readings)

            relative_log2_readings_by_protein[
                protein
            ] = relative_log2_normalised_readings

            results[protein][PROTEIN_ABUNDANCES][RAW] = readings
            results[protein][PROTEIN_ABUNDANCES][RAW][ABUNDANCE_AVERAGE] = raw_across_replicates_by_stage
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MEDIAN] = normalised_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MEDIAN][ABUNDANCE_AVERAGE] = normalised_means_across_replicates_by_stage
            # TODO - confirm the output later calculations are as they should be after this
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN] = copy.deepcopy(relative_log2_normalised_readings)
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN][ABUNDANCE_AVERAGE] = relative_log2_normalised_means_across_replicates_by_stage
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][ZERO_MAX] = min_max_normalised_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][ZERO_MAX][ABUNDANCE_AVERAGE] = min_max_normalised_means_across_replicates_by_stage
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MIN_MAX] = level_two_normalised_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MIN_MAX][ABUNDANCE_AVERAGE] = level_two_normalised_means_across_replicates_by_stage
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_PALBO] = arrest_log2_normalised_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_PALBO][ABUNDANCE_AVERAGE] = arrest_log2_normalised_means_across_replicates_by_stage
            results[protein][PROTEIN_ABUNDANCES][IMPUTED] = imputed_readings
            results[protein][PROTEIN_ABUNDANCES][IMPUTED][ABUNDANCE_AVERAGE] = imputed_across_replicates_by_stage
            results[protein][METRICS][LOG2_MEAN] = log2_mean_metrics
            results[protein][METRICS][ZERO_MAX] = zero_max_mean_metrics

        fisher_stats = self._calculate_fisher(relative_log2_readings_by_protein)

        # TODO - fisher stats vary wildly from ICR, investigate
        # print("++++++ FISHER STATS")
        # print(fisher_stats[FOCUS_PROTEIN])
        # return

        # TODO - rename this
        prot_anova_info: dict = {}
        for protein in anovas.keys():
            prot_anova_info[protein] = {}

            prot_anova_info[protein]["p_value"] = anovas[protein]["p_value"]

        # TODO - converting to a dataframe seems excessive. Find an alternative.
        prot_anova_info_df = pd.DataFrame(prot_anova_info).T
        prot_anova_info_df["q_value"] = self.p_adjust_bh(prot_anova_info_df["p_value"])

        prot_anova_info = prot_anova_info_df.to_dict("index")

        for protein in raw_readings:
            results[protein][METRICS][LOG2_MEAN][ANOVA] = anovas[protein]

            # ANOVA q values
            q_value = 1

            if protein in prot_anova_info:
                q_value = prot_anova_info[protein]["q_value"]

            results[protein][METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = q_value

            # Fisher
            fisher = {"G_statistic": 1, "p_value": 1, "frequency": 1, "q_value": 1}

            if protein in fisher_stats:
                fisher = fisher_stats[protein]

            results[protein][METRICS][LOG2_MEAN][FISHER_G] = fisher

        print(f"Number of proteins: {num_proteins}")
        print(json.dumps(results[FOCUS_PROTEIN]))
        return

    # TODO - lifted from ICR, rename
    def addPhosphoMetrics(self, time_course_phospho):
        """
        Add all metrics for each phosphosite.
        """
        for uniprot_accession in time_course_phospho:
            for site in time_course_phospho[uniprot_accession]["phosphorylation_abundances"]:
                norm_abundances = time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["position_abundances"]["normalised"]

                time_course_phospho[uniprot_accession]["phosphorylation_abundances"][
                    site
                ]["metrics"] = {
                    "log2_mean": self._calcAbundanceMetrics(
                        norm_abundances["log2_mean"], uniprot_accession),
                    "0-max": self._calcAbundanceMetrics(
                        norm_abundances["0-max"], uniprot_accession
                    )
                }
                # ANOVA
                p_value, f_statistic = self._calcANOVA(norm_abundances["log2_mean"])
                time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["ANOVA"] = {}
                time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["ANOVA"]["p_value"] = p_value
                time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["ANOVA"]["F_statistics"] = f_statistic
    
        # Fisher G Statistic
        time_course_fisher_dict = self._calcFisherG(time_course_phospho, "log2_mean", raw = False, phospho = True)

        # Corrected q values - Phospho
        # 1) Create a dataframe with the desired regression info
        phospho_anova_info = {}
        for uniprot_accession in time_course_phospho:
            for site in time_course_phospho[uniprot_accession]["phosphorylation_abundances"]:
                phospho_key = uniprot_accession + "_" + time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]['phosphorylation_site']
                if phospho_key not in phospho_anova_info:
                    phospho_anova_info[phospho_key] = {}
                phospho_anova_info[phospho_key]['p_value'] = time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["ANOVA"]["p_value"]
    
        phospho_anova_info_df = pd.DataFrame(phospho_anova_info)
        phospho_anova_info_df = phospho_anova_info_df.T
        # 2) Regression ANOVA q values
        phospho_anova_info_df['q_value'] = p_adjust_bh(phospho_anova_info_df['p_value'])
        # 3) Turn dataframe into a dictionary
        phospho_anova_info = phospho_anova_info_df.to_dict('index')
        # 4) Add Regression info in time_course_phospho dictionary
        for uniprot_accession in time_course_phospho:
            for site in time_course_phospho[uniprot_accession]["phosphorylation_abundances"]:
                site_key = uniprot_accession + "_" + time_course_phospho[uniprot_accession]['phosphorylation_abundances'][site]['phosphorylation_site']
                # ANOVA q values
                if site_key in phospho_anova_info:
                    time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["ANOVA"]["q_value"] = phospho_anova_info[site_key]['q_value']
                else:
                    time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["ANOVA"]["q_value"] = 1
                # # Fisher
                # if site_key in time_course_fisher_dict:
                #     time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["Fisher_G"] = time_course_fisher_dict[site_key]
                # else:
                #     time_course_phospho[uniprot_accession]["phosphorylation_abundances"][site]["metrics"]["log2_mean"]["Fisher_G"] = {'G_statistic': 1, 'p_value': 1, 'frequency': 1, 'q_value': 1}

        return time_course_phospho


    def _calculate_fisher(
        self,
        combined_time_course_info,
        phospho=False,
        phospho_ab=False,
        phospho_reg=False,
    ):
        # TODO - remove this comment
        # These calls are always the same in ICR
        # norm_method = "log2_mean"
        # raw = False

        if phospho:
            raise Exception("Phospho not implemented yet.")

        if phospho_ab:
            raise Exception("Phospho_ab not implemented yet.")

        if phospho_reg:
            raise Exception("Phospho_reg not implemented yet.")

        time_course_fisher = self.create_readings_dataframe(combined_time_course_info)

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

        time_course_fisher["G_statistic"] = g_stats
        time_course_fisher["p_value"] = p_values
        time_course_fisher["frequency"] = frequencies

        time_course_fisher["q_value"] = self.p_adjust_bh(time_course_fisher["p_value"])

        # Return only the Fisher columns
        fisher_cols = ["G_statistic", "p_value", "frequency", "q_value"]
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
    def create_readings_dataframe(self, readings):
        # TODO - raw is always False, ignore all raw code
        # TODO - norm_method is always log2_mean, ignore all norm_method code
        abundance_table = {}
        final_protein = None

        # TODO - maybe this and other similar loops could be changed to use items()?
        #   That way the value variable could be used in each successive loop, e.g.
        #   for protein, p_readings in readings.items():
        #       for replicate_name, rn_readings = p_readings.items():
        #           etc.
        for protein in readings:
            final_protein = protein

            abundance_table[protein] = {}

            for replicate_name in readings[protein]:
                for stage_name in readings[protein][replicate_name].keys():
                    rep_stage_name = f"{replicate_name}_{stage_name}"
                    abundance_table[protein][rep_stage_name] = readings[protein][
                        replicate_name
                    ][stage_name]

        time_course_abundance_df = pd.DataFrame(abundance_table)
        time_course_abundance_df = time_course_abundance_df.T

        new_cols = []

        replicate_names = list(readings[final_protein].keys())

        for stage_name in readings[final_protein][replicate_names[0]]:
            for rn in replicate_names:
                new_cols.append(f"{rn}_{stage_name}")

        time_course_abundance_df = time_course_abundance_df[new_cols]

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
    def _calcANOVA(self, readings: dict):
        """
        Groups by time point and performs ANOVA between all time point groups.
        """
        # Defaults if not enough replicate
        # TODO - why these values?
        p_value = 1
        f_statistic = 1

        first_replicate_name = list(readings.keys())[0]

        # TODO - change to stage_names
        timepoints_1 = []

        try:
            for stage_name in readings[first_replicate_name]:
                timepoints_1.append(self._tp(stage_name, readings))
    
            # TODO - is this necessary?
            timepoints = [x for x in timepoints_1 if x != []]

            # Protein info in at least 2 reps:
            # TODO - why is that comment above? There's no mention of two reps.
            if len(timepoints) != 0:
                one_way_anova = f_oneway(*timepoints)
                f_statistic = one_way_anova[0].item()
                p_value = one_way_anova[1].item()
                if np.isnan(p_value):
                    p_value = 1
        except Exception as e:
            print("++++ ERROR CALCULATING ANOVA")
            print(e)

        return {
            "p_value": p_value,
            "f_statistic": f_statistic
        }

    def _calculate_metrics(
        self,
        readings: dict,
        readings_averages: dict,
    ):
        metrics = {}

        abundances = []
        abundance_averages = readings_averages
        # TODO - how does ICR cope with Nones but not this? Does it use imputed values?
        abundance_averages_list = [
            val for val in abundance_averages.values() if val is not None
        ]
        # abundance_averages_list = list(abundance_averages.values())

        for replicate_name in readings:
            for stage_name in readings[replicate_name]:
                abundance = readings[replicate_name][stage_name]

                if abundance is not None:
                    # TODO - why does this add all reps for standard deviation calculation?
                    abundances.append(abundance)

        std_dev = None

        if len(abundances) > 1:
            std_dev = statistics.stdev(abundances)
        else:
            print("+++ NO ABUNDANCES FOR STANDARD DEVIATION")
            print(readings)

        try:
            _, residuals, *_ = np.polyfit(
                range(len(abundance_averages)),
                abundance_averages_list,
                2,
                full=True,
            )

            if len(residuals) == 0:
                # TODO - why 5?
                # TODO - 5 for ICR, but what about others?
                # eg Q9HBL0 {'G2_2': 0.4496, 'G2/M_1': 0.7425, 'M/Early G1': 1.0}
                residual = 5
            else:
                residual = residuals[0]

            # TODO - find out what this all means
            r_squared = self._polyfit(
                range(0, len(abundance_averages)), abundance_averages_list, 2
            )
            # TODO - find out what this all means
            max_fold_change = max(abundance_averages_list) - min(
                abundance_averages_list
            )

            # TODO - what does all this mean?
            base_metrics = {
                # Why is variance average 2 sig fig in ICR but several here?
                "variance": moment(abundance_averages_list, moment=2),
                "skewness": moment(abundance_averages_list, moment=3),
                "kurtosis": moment(abundance_averages_list, moment=4),
                "peak": max(abundance_averages, key=abundance_averages.get),
                "max_fold_change": max_fold_change,
                "residuals": residual,
                "R_squared": r_squared,
            }

            # If we have info for the protein in at least 2 replicates
            if len(readings) >= 2:
                curve_fold_change, curve_peak = self._calcCurveFoldChange(
                    # readings, protein.accession_number
                    readings
                )
                residuals_all, r_squared_all = self._calcResidualsR2All(readings)

                metrics = {
                    "standard_deviation": std_dev,
                    **{f"{k}_average": v for k, v in base_metrics.items()},
                    "residuals_all": residuals_all,
                    "R_squared_all": r_squared_all,
                    "curve_fold_change": curve_fold_change,
                    "curve_peak": curve_peak,
                }

        except Exception as e:
            print("Exception in _calculate_metrics")
            print(e)

        return metrics

    # TODO - this is straight up lifted from ICR, change and ideally use a library
    def _calcResidualsR2All(self, readings):
        """
        Calculate the residuals and the R squared for all the abundances from all the replicates for a protein.
        """
        residuals_all = None
        r_squared_all = None

        x, y, _ = self._generate_xs_ys(readings)

        if len(x) == len(y):
            p = np.poly1d(np.polyfit(x, y, 2))
            curve_abundances = p(x)
            residuals_all = np.polyfit(x, y, 2, full=True)[1][0]
            residuals_all = residuals_all.item()
            r_squared_all = round(r2_score(y, curve_abundances), 2)

        return residuals_all, r_squared_all

    # TODO - this is straight up lifted from ICR, change and ideally use a library
    def _calcCurveFoldChange(self, readings):
        """
        Calculates the curve_fold_change and curve peaks for the three or two replicates normalised abundance for each protein.
        """
        curve_fold_change = None
        curve_peak = None

        x, y, stage_names_map = self._generate_xs_ys(readings)

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

    # TODO - why does ICR not need this?
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
                f"Formatting for {protein_no}, {protein_reading.protein.accession_number}",
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
                f"Formatting for {mod_no}, protein {phospho_reading.phospho.protein.accession_number} mod {phospho_reading.phospho.mod}",
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

    def _calculate_means_across_replicates_by_stage(
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
        readings: dict,
    ):
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
                    "Can't create median with no abundances for protein {protein_reading.protein.name}"
                )

            median = statistics.median(readings)

            stage_name_medians[column_name.sample_stage.name] = median

        return stage_name_medians

    def _count_logger(self, i: int, step: int, output: str):
        if i % step == 0:
            logger.info(output)

    def _round(self, value):
        return round(value, 4)

    def _generate_xs_ys(self, readings):
        final_replicate_name = list(readings.keys())[-1]

        stage_names_map = {}

        for i, stage_name in enumerate(readings[final_replicate_name].keys()):
            stage_names_map[stage_name] = i

        x = []
        for stage_name in readings.get(final_replicate_name, {}):
            x.append(stage_names_map[stage_name])
        x.sort()

        y = []
        for stage_name in stage_names_map:
            value = readings.get(final_replicate_name, {}).get(stage_name)
            if value is not None:
                y.append(value)

        return x, y, stage_names_map