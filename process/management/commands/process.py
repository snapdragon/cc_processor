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

from process.models import ColumnName, Project, Phospho, ProteinReading, Replicate, PhosphoReading

# TODO - make this configurable by flag?
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

# TODO - move constants elsewhere
# FOCUS_PROTEIN_ACCESSION_NUMBER = "Q09666"
FOCUS_PROTEIN_ACCESSION_NUMBER = "O95714"
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
P_VALUE = "p_value"
F_STATISTICS = "F_statistics"
PROTEIN_OSCILLATION_ABUNDANCES = "protein_oscillation_abundances"
ABUNDANCE_AVERAGE = "abundance_average"

POSITION_ABUNDANCES = "position_abundances"
FOCUS_MOD = "2928"
PHOSPHORYLATION_ABUNDANCES = "phosphorylation_abundances"
PHOSPHO_REGRESSION = "phospho_regression"
PHOSPHORYLATION_SITE = "phosphorylation_site"
Q_VALUE = Q_VALUE

TOTAL_PROTEIN_INDEX_FILE = "total_protein_index.json"

# TOOD - what is this file? Where did it come from?
# TODO - whatever it is, a lot of it could be trimmed, and maybe put in the DB
with open(f"./data/{TOTAL_PROTEIN_INDEX_FILE}") as outfile:
    index_protein_names = json.load(outfile)


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

        results = self._proteo(project, replicates, protein_readings, column_names, with_bugs)

        # phospho_readings = PhosphoReading.objects.filter(phospho__protein__project=project)
        phospho_readings = PhosphoReading.objects.filter(phospho__protein__project=project)[:1000]
        phosphos = Phospho.objects.filter(protein__project=project)

        phospho_results = self._phospho(project, replicates, phospho_readings, column_names, phosphos, with_bugs)

        # print("++++ NO PROTEINS BEFORE")
        # print(len(list(results.keys())))

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

        # print("++++ NO PROTEINS")
        # print(len(list(results.keys())))

        self._addProteinOscillations(results, with_bugs)



    def _phospho(self, project, replicates, phospho_readings, column_names, phosphos, with_bugs: bool):
        logger.info("Processing phosphoproteome")

        # TODO - convert this in to a DB query to save having to run it manually?
        readings_by_rep_stage = self._format_phospho_readings(phospho_readings)
        phosphosites = self._get_phosphosites(phosphos)

        # raw_readings = self._qc_phospho_readings(readings_by_rep_stage)

        del(phospho_readings)
        del(phosphos)

        raw_readings = readings_by_rep_stage

        medians = self._calculate_phospho_medians(raw_readings)

        num_proteins = 0

        relative_log2_readings_by_protein = {}
        anovas = {}
        results = {}

        # TODO - make this a flag
        # TODO - this is a duplicate
        # FOCUS_PROTEIN = Protein.objects.get(
        #     accession_number=FOCUS_PROTEIN_ACCESSION_NUMBER, project__name=project.name
        # )

        for protein in raw_readings.keys():
            results[protein] = {}

            for mod, readings in raw_readings[protein].items():
                num_proteins += 1

                self._count_logger(
                    num_proteins,
                    10000,
                    f"Formatting for {num_proteins}, {protein.accession_number}",
                )

                results[protein][mod] = {
                    PHOSPHORYLATION_SITE: phosphosites[mod],
                    POSITION_ABUNDANCES: {
                        RAW: {},
                        NORMALISED: {},
                        IMPUTED: {}
                    },
                    METRICS: {}
                }

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
                    normalised_readings
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

                # TODO - why is this called median when it's from normalised?
                median_averages = (
                    self._calculate_means(
                        normalised_readings, with_bugs
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
                )

                zero_max_mean_metrics = self._calculate_metrics(
                    zero_max_readings,
                    zero_max_averages,
                )

                anovas[protein] = self._calcANOVA(log2_readings)

                relative_log2_readings_by_protein[
                    protein
                ] = log2_readings

                # TODO - these calculations and these results can be combined in one function
                #   with proteo
                results[protein][mod][POSITION_ABUNDANCES][RAW] = readings
                results[protein][mod][POSITION_ABUNDANCES][RAW][ABUNDANCE_AVERAGE] = raw_averages
                # TODO - confirm the output later calculations are as they should be after this
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][LOG2_MEAN] = copy.deepcopy(log2_readings)
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][LOG2_MEAN][ABUNDANCE_AVERAGE] = log2_averages
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][MIN_MAX] = min_max_readings
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][MIN_MAX][ABUNDANCE_AVERAGE] = min_max_averages
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][MEDIAN] = normalised_readings
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][MEDIAN][ABUNDANCE_AVERAGE] = normalised_averages
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][ZERO_MAX] = zero_max_readings
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][ZERO_MAX][ABUNDANCE_AVERAGE] = zero_max_averages
                # TODO - this should be called 'arrest'
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][LOG2_PALBO] = arrest_readings
                results[protein][mod][POSITION_ABUNDANCES][NORMALISED][LOG2_PALBO][ABUNDANCE_AVERAGE] = arrest_averages
                results[protein][mod][POSITION_ABUNDANCES][IMPUTED] = imputed_readings
                results[protein][mod][POSITION_ABUNDANCES][IMPUTED][ABUNDANCE_AVERAGE] = imputed_averages
                results[protein][mod][METRICS][LOG2_MEAN] = log2_mean_metrics
                results[protein][mod][METRICS][ZERO_MAX] = zero_max_mean_metrics

        # # Kinase Consensus Prediction
        # # TODO - lifted, change
        # # TODO - what does all this do?
        # # TODO - put back in and fix
        # for protein in raw_readings:
        #     for mod in raw_readings[protein]:
        #         # only for certain phosphorylation sites
        #         # TODO - why? What does the lack of dash mean?
        #         if mod.find("-") == -1:
        #             print(f"+++++ NO DASH {protein} {mod}")
        #             # Add PepTools Phospho annotations
        #             mod_result = results[protein][mod]

        #             if 'phospho' in index_protein_names[protein.accession_number] and mod in index_protein_names[protein.accession_number]['phospho']:
        #                 mod_result['PepTools_annotations'] = index_protein_names[protein.accession_number]['phospho'][mod]

        #             mod_result['kinase_prediction'] = {}

        #             phospho_kinases_class = self._getConsensusKinasePred(protein.accession_number, mod_result)
        #             mod_result['kinase_prediction']['peptide_seq'] = phospho_kinases_class['peptide_seq']
        #             mod_result['kinase_prediction']['consenus_motif_match'] = phospho_kinases_class['kinase_motif_match']

        # print("++++++ RESULTS")
        # print(json.dumps(results[FOCUS_PROTEIN][FOCUS_MOD]))
        # # print(results.keys())
        # exit()

        return results


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
        # FOCUS_PROTEIN = Protein.objects.get(
        #     accession_number=FOCUS_PROTEIN_ACCESSION_NUMBER, project__name=project.name
        # )

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
        #     self._calculate_means(
        #         protein_readings, with_bugs
        #     )
        # )

        # N.B. protein_readings_by_rep_stage is not the same structure as protein_readings.
        #   protein_readings is just a list of ProteinReading objects. normalised_protein_readings
        #   is a dict with Protein object keys. The values is a dict of replicate name keys
        #   with a dict of sample stage names and abundances as keys.
        # TODO - convert this in to a DB query to save having to run it manually?
        readings_by_rep_stage = self._format_protein_readings(protein_readings)

        del(protein_readings)

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
                PHOSPHORYLATION_ABUNDANCES: {},
                METRICS: {}
            }

            # TODO - check whether all means calculations need with-bugs
            # TODO - change all these to 'raw_averages' and the like
            raw_averages = (
                self._calculate_means(
                    readings, with_bugs
                )
            )

            # firstLevelNormalisationProteomics
            normalised_readings = self._calculate_first_level_normalisation(
                readings, medians
            )

            normalised_averages = (
                self._calculate_means(
                    normalised_readings, with_bugs
                )
            )

            # calclog2PalboNormalisation
            arrest_readings = self._calculate_arrest_log2_normalisation(
                normalised_readings, project
            )

            arrest_averages = (
                self._calculate_means(
                    arrest_readings, with_bugs
                )
            )

            # calclog2RelativeAbundance
            log2_readings = (
                self._calculate_relative_log2_normalisation(normalised_readings)
            )

            log2_averages = (
                self._calculate_means(
                    log2_readings, with_bugs
                )
            )

            # normaliseData
            min_max_readings = self._calculate_level_two_normalisation(
                log2_readings
            )

            min_max_averages = (
                self._calculate_means(
                    min_max_readings, with_bugs
                )
            )

            zero_max_readings = self._calculate_level_two_normalisation(
                normalised_readings, True
            )

            zero_max_averages = (
                self._calculate_means(
                    zero_max_readings, with_bugs
                )
            )

            imputed_readings = self._impute(
                min_max_readings, replicates, column_names
            )

            imputed_averages = (
                self._calculate_means(
                    imputed_readings, with_bugs, imputed=True
                )
            )

            log2_mean_metrics = self._calculate_metrics(
                log2_readings,
                log2_averages,
            )

            zero_max_mean_metrics = self._calculate_metrics(
                zero_max_readings,
                zero_max_averages,
            )

            anovas[protein] = self._calcANOVA(log2_readings)

            relative_log2_readings_by_protein[
                protein
            ] = log2_readings

            results[protein][PROTEIN_ABUNDANCES][RAW] = readings
            results[protein][PROTEIN_ABUNDANCES][RAW][ABUNDANCE_AVERAGE] = raw_averages
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MEDIAN] = normalised_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MEDIAN][ABUNDANCE_AVERAGE] = normalised_averages
            # TODO - confirm the output later calculations are as they should be after this
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN] = copy.deepcopy(log2_readings)
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN][ABUNDANCE_AVERAGE] = log2_averages
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][ZERO_MAX] = zero_max_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][ZERO_MAX][ABUNDANCE_AVERAGE] = zero_max_averages
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MIN_MAX] = min_max_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][MIN_MAX][ABUNDANCE_AVERAGE] = min_max_averages
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_PALBO] = arrest_readings
            results[protein][PROTEIN_ABUNDANCES][NORMALISED][LOG2_PALBO][ABUNDANCE_AVERAGE] = arrest_averages
            results[protein][PROTEIN_ABUNDANCES][IMPUTED] = imputed_readings
            results[protein][PROTEIN_ABUNDANCES][IMPUTED][ABUNDANCE_AVERAGE] = imputed_averages
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

            prot_anova_info[protein][P_VALUE] = anovas[protein][P_VALUE]

        # TODO - converting to a dataframe seems excessive. Find an alternative.
        prot_anova_info_df = pd.DataFrame(prot_anova_info).T
        prot_anova_info_df[Q_VALUE] = self.p_adjust_bh(prot_anova_info_df[P_VALUE])

        prot_anova_info = prot_anova_info_df.to_dict("index")

        for protein in raw_readings:
            results[protein][METRICS][LOG2_MEAN][ANOVA] = anovas[protein]

            # ANOVA q values
            q_value = 1

            if protein in prot_anova_info:
                q_value = prot_anova_info[protein][Q_VALUE]

            results[protein][METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = q_value

            # Fisher
            fisher = {"G_statistic": 1, P_VALUE: 1, "frequency": 1, Q_VALUE: 1}

            if protein in fisher_stats:
                fisher = fisher_stats[protein]

            results[protein][METRICS][LOG2_MEAN][FISHER_G] = fisher

        # print(f"Number of proteins: {num_proteins}")
        # print(json.dumps(results[FOCUS_PROTEIN]))
        return results

    # TODO - lifted from ICR, rename variables
    # TODO - check all comments are not from ICR
    def _calculate_phospho_metrics(self, readings, log2_readings, log2_averages, zero_max_readings, zero_max_averages):
        """
        Add all metrics for each phosphosite.
        """
        metrics = {}

        # TODO - finish all these
        # TODO - looks like some duplication, simplify    
        # # Fisher G Statistic
        # time_course_fisher_dict = self._calculate_fisher(time_course_phospho, phospho = True)

        # Corrected q values - Phospho
        # 1) Create a dataframe with the desired regression info
        # phospho_anova_info = {}
        # for uniprot_accession in time_course_phospho:
        #     for site in time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES]:
        #         phospho_key = uniprot_accession + "_" + time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]
        #         if phospho_key not in phospho_anova_info:
        #             phospho_anova_info[phospho_key] = {}
        #         phospho_anova_info[phospho_key]['p_value'] = time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][METRICS][LOG2_MEAN][ANOVA][P_VALUE]
    
        # phospho_anova_info_df = pd.DataFrame(phospho_anova_info)
        # phospho_anova_info_df = phospho_anova_info_df.T
        # # 2) Regression ANOVA q values
        # phospho_anova_info_df['q_value'] = p_adjust_bh(phospho_anova_info_df['p_value'])
        # # 3) Turn dataframe into a dictionary
        # phospho_anova_info = phospho_anova_info_df.to_dict('index')
        # # 4) Add Regression info in time_course_phospho dictionary
        # for uniprot_accession in time_course_phospho:
        #     for site in time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES]:
        #         site_key = uniprot_accession + "_" + time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]
        #         # ANOVA q values
        #         if site_key in phospho_anova_info:
        #             time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = phospho_anova_info[site_key]['q_value']
        #         else:
        #             time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][METRICS][LOG2_MEAN][ANOVA][Q_VALUE] = 1
        #         # # Fisher
        #         # if site_key in time_course_fisher_dict:
        #         #     time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][METRICS][LOG2_MEAN]["Fisher_G"] = time_course_fisher_dict[site_key]
        #         # else:
        #         #     time_course_phospho[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][METRICS][LOG2_MEAN]["Fisher_G"] = {'G_statistic': 1, 'p_value': 1, 'frequency': 1, 'q_value': 1}

        return metrics


    def _calculate_fisher(
        self,
        results,
        phospho=False,
        phospho_ab=False,
        phospho_reg=False,
    ):
        # TODO - remove this comment
        # These calls are always the same in ICR
        # norm_method = LOG2_MEAN
        # raw = False

        time_course_fisher = self._create_results_dataframe(results, phospho, phospho_ab, phospho_reg)

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
        time_course_fisher[P_VALUE] = p_values
        time_course_fisher["frequency"] = frequencies

        time_course_fisher[Q_VALUE] = self.p_adjust_bh(time_course_fisher[P_VALUE])

        # Return only the Fisher columns
        fisher_cols = ["G_statistic", P_VALUE, "frequency", Q_VALUE]
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
    def _create_results_dataframe(self, results, phospho, phospho_ab, phospho_reg):
        # TODO - raw is always False, ignore all raw code
        # TODO - norm_method is always log2_mean, ignore all norm_method code
        abundance_table = {}
        final_protein = None
        final_mod = None

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
                # print("++++ PRS")
                # print(protein)
                # print(prs[PHOSPHORYLATION_ABUNDANCES])

                for mod in prs[PHOSPHORYLATION_ABUNDANCES]:
                    final_mod = mod

                    protein_abundances_all = prs[PHOSPHORYLATION_ABUNDANCES][mod][POSITION_ABUNDANCES]

                    # print("++++ PAA")
                    # print(protein_abundances_all)
                    # exit()

                    mod_key = protein.accession_number + "_" + prs[PHOSPHORYLATION_ABUNDANCES][mod][PHOSPHORYLATION_SITE]

                    # print("+++++ LEN")
                    # print(len(protein_abundances_all))

                    if len(protein_abundances_all) != 0:
                        # TODO - tidy this
                        protein_abundances = protein_abundances_all[NORMALISED][LOG2_MEAN]

                        # print("++++++++ FOO")
                        # # print(protein_abundances)

                        # print(prs[PHOSPHORYLATION_ABUNDANCES][mod].keys())

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

                        # TODO - make generic
                        # TODO - tidy this
                        for abundance_rep in ['abundance_One','abundance_Two']:
                            if mod_key not in abundance_table:
                                abundance_table[mod_key] = {}

                            # TODO - not in original, why here?
                            if not protein_abundances.get(abundance_rep):
                                continue

                            for timepoint in protein_abundances[abundance_rep]:
                                rep = "_".join(abundance_rep.split("_", 1)[1].split("_")[:2])
                                rep_timepoint = rep + "_" + timepoint
                                abundance_table[mod_key][rep_timepoint] = protein_abundances[abundance_rep][timepoint]

            else:
                abundance_table[protein] = {}

                for replicate_name in results[protein]:
                    for stage_name in results[protein][replicate_name].keys():
                        rep_stage_name = f"{replicate_name}_{stage_name}"
                        abundance_table[protein][rep_stage_name] = results[protein][
                            replicate_name
                        ][stage_name]

        # if phospho:
        #     print("++++ A TABLE")
        #     print(abundance_table)
        #     exit()

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

        # print("++++++ R NAMES")
        # print(replicate_names)

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

        # print("++++ DF")
        # print(time_course_abundance_df)
        # exit()

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

        # TODO - don't attempt if less than two values, it errors
        try:
            for stage_name in readings[first_replicate_name]:
                timepoints_1.append(self._tp(stage_name, readings))
    
            # TODO - is this necessary?
            timepoints = [x for x in timepoints_1 if x != []]

            # f_oneway needs at least 2 stage names
            if len(timepoints) > 1:
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
            print("NO ABUNDANCES FOR STANDARD DEVIATION")
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
                f"Formatting for {phospho.mod}",
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
                    "Can't create median with no abundances for protein {protein_reading.protein.accession_number}"
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
    def _addProteinOscillations(self, results, with_bugs):
        # TODO - check all logging statements for similarity to ICR
        logger.info("Adding Protein Oscillation Normalised Abundances")

        for protein in results:
            # If we have info both in protein and in phospho level
            prpa = results[protein][PHOSPHORYLATION_ABUNDANCES]

            # print("+++++ BAR")
            # print(len(prpa))
            # print(prpa)
            # print(len(results[protein][PROTEIN_ABUNDANCES][RAW]))
            # print(results[protein][PROTEIN_ABUNDANCES][RAW])

            if len(prpa) != 0 and len(results[protein][PROTEIN_ABUNDANCES][RAW]) != 0:
                # print("+++++ ADDING")
                # exit()

                # Add the Protein Oscillation Normalised Abundances
                for mod in prpa:
                    # TODO - could be initialised in loop below
                    prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES] = {ZERO_MAX:{}, LOG2_MEAN:{}}
                    prpampoa = prpa[mod][PROTEIN_OSCILLATION_ABUNDANCES]

                    for norm_method in [ZERO_MAX, LOG2_MEAN]:
                        protein_oscillation_abundances = prpampoa[norm_method]

                        phospho_oscillations = self._calculateProteinOscillationAbundances(
                            results, norm_method, protein, mod
                        )
                        # TODO - make generic
                        for rep in ["One", "Two"]:
                            if rep in phospho_oscillations:
                                key = "abundance_" + rep
                                protein_oscillation_abundances[key] = phospho_oscillations[rep]

                        protein_oscillation_abundances[ABUNDANCE_AVERAGE] = self._calculate_means(
                            protein_oscillation_abundances, imputed=False, with_bugs=with_bugs)

                    # if we have info in Protein Oscillation Normalised Abundances
                    if (
                        len(prpampoa[LOG2_MEAN]) > 1
                    ):
                        for norm_method in [ZERO_MAX, LOG2_MEAN]:
                            protein_oscillation_abundances = prpampoa[norm_method]

                            # Metrics
                            prpampoa[norm_method][METRICS] = self._calcAbundanceMetrics(
                                protein_oscillation_abundances, protein)

                        # ANOVA
                        norm_abundances = prpampoa[LOG2_MEAN]

                        anovas = self._calcANOVA(norm_abundances)

                        prpampoa[LOG2_MEAN][METRICS][ANOVA] = anovas

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

        prot_phospho_info_df = pd.DataFrame(prot_phospho_info)
        prot_phospho_info_df = prot_phospho_info_df.T
        # 2) Protein-Phospho ANOVA q values
        prot_phospho_info_df['q_value'] = self.p_adjust_bh(prot_phospho_info_df['p_value'])
        # 3) Turn dataframe into a dictionary
        prot_phospho_info = prot_phospho_info_df.to_dict('index')

        # 4) Add Protein-Phospho info in combined_time_course_info dictionary
        # TODO - tidy up, two loops not necessary?
        for protein in results:
            if len(results[protein][PHOSPHORYLATION_ABUNDANCES]) != 0 and len(results[protein][PROTEIN_ABUNDANCES][RAW]) != 0:
                for site in results[protein][PHOSPHORYLATION_ABUNDANCES]:
                    if PROTEIN_OSCILLATION_ABUNDANCES in results[protein][PHOSPHORYLATION_ABUNDANCES][site]:
                        site_key = protein.accession_number + "_" + results[protein][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]
                        # ANOVA q values
                        # TODO - tidy
                        if site_key in prot_phospho_info:
                            results[protein][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = prot_phospho_info[site_key]['q_value']
                        else:
                            results[protein][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = 1
                        # # Fisher
                        # if site_key in time_course_fisher_dict:
                        #     results[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS]["Fisher_G"] = time_course_fisher_dict[site_key]
                        # else:
                        #     results[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PROTEIN_OSCILLATION_ABUNDANCES][LOG2_MEAN][METRICS]["Fisher_G"] = {'G_statistic': 1, 'p_value': 1, 'frequency': 1, 'q_value': 1}

        return results






# def addPhosphoRegression(combined_time_course_info):
#     """
#     Normalise the phospho abundance on the protein abundance
#     Calculates and adds all the Regressed Phospho Abundance and their metrics for each phosphosite.
#     Phospho = Dependent = Y
#     Protein = Independent = X
#     Linear Model => y = ax + b
#     Residuals = Y - Y_predict 
#     """
#     logger.info("Adding Phospho Normalised on Protein Abundances - Regression")

#     for uniprot_accession in combined_time_course_info:
#         # If we have info both in protein and in phospho level
#         if len(combined_time_course_info[uniprot_accession][
#                 PHOSPHORYLATION_ABUNDANCES]) != 0 and len(combined_time_course_info[uniprot_accession][
#                 PROTEIN_ABUNDANCES][RAW]) != 0:
            
#             for phosphosite in combined_time_course_info[uniprot_accession][
#                 PHOSPHORYLATION_ABUNDANCES
#             ]:  
#                 rep1_phosho = combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                     phosphosite][POSITION_ABUNDANCES][NORMALISED][LOG2_MEAN]['abundance_rep_1']
#                 rep2_phosho = combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                     phosphosite][POSITION_ABUNDANCES][NORMALISED][LOG2_MEAN]['abundance_rep_2']
#                 rep1_prot = combined_time_course_info[uniprot_accession][PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN]['abundance_rep_1']
#                 rep2_prot = combined_time_course_info[uniprot_accession][PROTEIN_ABUNDANCES][NORMALISED][LOG2_MEAN]['abundance_rep_2']
                
#                 # If we don't have missing values in the protein and phospho abundances
#                 if len(rep1_phosho) != 8 or len(rep1_prot) != 8 or len(rep2_phosho) != 8 or len(rep2_prot) != 8:
#                     continue
                
#                 # Add the Regressed Phospho Normalised Abundances
#                 combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                     phosphosite][PHOSPHO_REGRESSION] = {ZERO_MAX:{}, LOG2_MEAN:{}}

#                 # converting dictionary values to list for both replicates
#                 phospho_Y_list = list(rep1_phosho.values()) + list(rep2_phosho.values())
#                 protein_X_list = list(rep1_prot.values()) + list(rep2_prot.values())
#                 # converting list to array
#                 prot_x = np.asarray(protein_X_list).reshape(len(protein_X_list), 1)
#                 phospho_y = np.asarray(phospho_Y_list).reshape(len(phospho_Y_list), 1)
#                 # Fit the linear model
#                 model = linear_model.LinearRegression().fit(prot_x,phospho_y)
#                 # Predict new Phospho Values
#                 y_pred = model.predict(prot_x)
#                 # Calculate Residuals
#                 residuals = (phospho_y - y_pred)
#                 # Create new regressed phospho abundances dictionaries
#                 replicates = ["abundance_rep_1", "abundance_rep_2"]
#                 res_dic = {}   
#                 for replicate in replicates:
#                     res_dic[replicate] = {}
#                     if replicate == "abundance_rep_1":
#                         for index, value in enumerate(residuals[0:8]):
#                             key = time_points[index]
#                             res_dic[replicate][key] = value[0]
#                     if replicate == "abundance_rep_2":
#                         for index, value in enumerate(residuals[8::]):
#                             key = time_points[index]
#                             res_dic[replicate][key] = value[0]

#                 phospho_regression = combined_time_course_info[
#                     uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite][
#                     PHOSPHO_REGRESSION]
                
#                 phospho_regression[LOG2_MEAN] = res_dic

#                 phospho_regression[LOG2_MEAN][
#                         ABUNDANCE_AVERAGE
#                     ] = calculateAverageRepAbundance(
#                         phospho_regression[LOG2_MEAN],
#                         imputed=False,
#                         phospho_oscillation=False,
#                     )

#                 # Add metrics
#                 # Calculate the protein - phospho vector correlation
#                 phospho_regression[ZERO_MAX][METRICS] = {}
#                 phospho_regression[ZERO_MAX][METRICS]["protein-phosho-correlation"] = stats.pearsonr(protein_X_list, phospho_Y_list)[0] # [0] to get the correlation coefficient, [1] = p-value
#                 # curve fold change phosphorylation/curve fold change protein for  0-max
#                 phospho_regression[ZERO_MAX][METRICS]["phosho-protein-cfc_ratio"] = combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                     phosphosite][METRICS]['0-max']['curve_fold_change'] / combined_time_course_info[uniprot_accession][METRICS]['0-max']['curve_fold_change']

#                 # ANOVA
#                 p_value, f_statistic = calcANOVA(phospho_regression[LOG2_MEAN])
#                 phospho_regression[LOG2_MEAN][METRICS] = {}
#                 phospho_regression[LOG2_MEAN][METRICS][ANOVA] = {}
#                 combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                     phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][P_VALUE] = p_value
#                 combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                     phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][F_STATISTICS] = f_statistic

#     # Fisher G Statistic - Phospho
#     time_course_fisher_dict = self._calculate_fisher(combined_time_course_info, phospho = True, phospho_ab = False, phospho_reg = True)
#     # Corrected q values - Phospho
#     # 1) Create a dataframe with the desired regression info
#     regression_info = {}
#     for uniprot_accession in combined_time_course_info:
#         if len(combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES]) != 0:
#             for site in combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES]:
#                 if 'phospho_regression' in combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site]:
#                     phospho_key = uniprot_accession + "_" + combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site][PHOSPHORYLATION_SITE]
#                     if phospho_key not in regression_info:
#                         regression_info[phospho_key] = {}
#                     regression_info[phospho_key]['p_value'] = combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][site]['phospho_regression']['log2_mean'][METRICS]['ANOVA']['p_value']

#     regression_info_df = pd.DataFrame(regression_info)
#     regression_info_df = regression_info_df.T
#     # 2) Regression ANOVA q values
#     regression_info_df['q_value'] = p_adjust_bh(regression_info_df['p_value'])
#     # 3) Turn dataframe into a dictionary
#     regression_info = regression_info_df.to_dict('index')
#     # 4) Add Regression info in combined_time_course_info dictionary
#     for uniprot_accession in combined_time_course_info:
#         if len(combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES]) != 0 and len(combined_time_course_info[uniprot_accession][PROTEIN_ABUNDANCES][RAW]) != 0:
#             for phosphosite in combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES]:
#                 if PHOSPHO_REGRESSION in combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite]:
#                     site_key = uniprot_accession + "_" + combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHORYLATION_SITE]
#                     # ANOVA q values
#                     if site_key in regression_info:
#                         combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = regression_info[site_key]['q_value']
#                     else:
#                         combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS][ANOVA][Q_VALUE] = 1
#                     # Fisher
#                     phospho_regression[LOG2_MEAN][METRICS]["Fisher_G"] = {}
#                     if site_key in time_course_fisher_dict:
#                         combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS]["Fisher_G"] = time_course_fisher_dict[site_key]
#                     else:
#                         combined_time_course_info[uniprot_accession][PHOSPHORYLATION_ABUNDANCES][
#                                     phosphosite][PHOSPHO_REGRESSION][LOG2_MEAN][METRICS]["Fisher_G"] = {'G_statistic': 1, 'p_value': 1, 'frequency': 1, 'q_value': 1}

#     return combined_time_course_info




    # TODO - lifted from ICR
    def _calculateProteinOscillationAbundances(self, combined_time_course_info, norm_method, uniprot_accession, phosphosite):
        """
        # Protein oscillation normalised abundances to the phospho part of the JSON too.
        # It would just be the normalised phopho abundances minus the normalised protein abundances per replicate.
        # And then the average of them.
        """
        phospho_oscillations = {}

        # protein exists
        if "One" in combined_time_course_info[uniprot_accession]["protein_abundances"]["normalised"][norm_method]:
            # site exists
            if ("One" in combined_time_course_info[uniprot_accession]["phosphorylation_abundances"][
                    phosphosite]["position_abundances"]["normalised"][norm_method]):

                protein_oscillation_abundances_rep_1 = {}
                protein_normed_abundance_rep_1 = combined_time_course_info[uniprot_accession][
                    "protein_abundances"]["normalised"][norm_method]["One"]
            
                phosho_normed_abundance_rep_1 = combined_time_course_info[uniprot_accession][
                    "phosphorylation_abundances"][phosphosite]["position_abundances"]["normalised"][norm_method][
                    "One"]
            
                for time_point in protein_normed_abundance_rep_1:
                    if time_point in phosho_normed_abundance_rep_1:
                        protein_oscillation_abundances_rep_1[time_point] = (
                            phosho_normed_abundance_rep_1[time_point]
                            - protein_normed_abundance_rep_1[time_point]
                        )
                phospho_oscillations["One"] = protein_oscillation_abundances_rep_1

        # Replicate 2
        if "Two" in combined_time_course_info[uniprot_accession]["protein_abundances"]["normalised"][norm_method]:
            if ("Two" in combined_time_course_info[uniprot_accession]["phosphorylation_abundances"][
                    phosphosite]["position_abundances"]["normalised"][norm_method]):
            
                protein_oscillation_abundances_rep_2 = {}
                protein_normed_abundance_rep_2 = combined_time_course_info[uniprot_accession][
                    "protein_abundances"]["normalised"][norm_method]["Two"]
                phosho_normed_abundance_rep_2 = combined_time_course_info[uniprot_accession][
                    "phosphorylation_abundances"][phosphosite]["position_abundances"]["normalised"][norm_method][
                    "Two"]
                for time_point in protein_normed_abundance_rep_2:
                    if time_point in phosho_normed_abundance_rep_2:
                        protein_oscillation_abundances_rep_2[time_point] = (
                            phosho_normed_abundance_rep_2[time_point]
                            - protein_normed_abundance_rep_2[time_point]
                            )
                phospho_oscillations["Two"] = protein_oscillation_abundances_rep_2
    
        return phospho_oscillations



    # TODO - lifted from ICR
    # TODO - there is duplicate here with another function
    # TODO - study this
    def _calcAbundanceMetrics(self, norm_abundances, uniprot_accession):
        """
        Calculates the moments and peaks for the average for the three replicates normalised abundance for each protein.
        """
        metrics = {}
        norm_abundance_average = norm_abundances["abundance_average"]
        norm_abundance_average_list = list(norm_abundance_average.values())

        norm_method = '0-max'

        # TODO - make generic
        time_points = [
	        "Palbo",
	        "Late G1_1",
	        "G1/S",
	        "S",
    	    "S/G2",
	        "G2_2",
    	    "G2/M_1",
	        "M/Early G1",
        ]

        abundance = []
        for timepoint in time_points:
            for rep in norm_abundances:
                if rep != "abundance_average":
                    if timepoint in norm_abundances[rep]:
                        abundance.append(norm_abundances[rep][timepoint])
    
        std = statistics.stdev(abundance)

        try:
            residuals_array = np.polyfit(
                range(0, len(norm_abundance_average)),
                norm_abundance_average_list,
                2,
                full=True,
            )[1]

            if len(residuals_array) == 0:
                # eg Q9HBL0 {'G2_2': 0.4496, 'G2/M_1': 0.7425, 'M/Early G1': 1.0}
                residuals = 5
            else:
                residuals = np.polyfit(
                range(0, len(norm_abundance_average)),
                norm_abundance_average_list,
                2,
                full=True,)[1][0]

            r_squared = np.polyfit(
                range(0, len(norm_abundance_average)), norm_abundance_average_list, 2
            )
            max_fold_change = max(norm_abundance_average.values()) - min(
                norm_abundance_average.values()
            )
            metrics = {
                "variance": moment(norm_abundance_average_list, moment=2),
                "skewness": moment(norm_abundance_average_list, moment=3),
                "kurtosis": moment(norm_abundance_average_list, moment=4),
                "peak": max(norm_abundance_average, key=norm_abundance_average.get),
                "max_fold_change": max_fold_change,
                "residuals": residuals,
                "R_squared": r_squared,
            }
            # if we have info for the protein in at least 2 replicates
            if len(norm_abundances) == 3:
                curve_fold_change, curve_peak = self._calcCurveFoldChange(
                    # norm_abundances, uniprot_accession
                    norm_abundances
                )
                residuals_all, r_squared_all = self._calcResidualsR2All(norm_abundances)
                metrics = {
                    "standard_deviation": std, 
                    "variance_average": round(moment(norm_abundance_average_list, moment=2),2),
                    "skewness_average": moment(norm_abundance_average_list, moment=3),
                    "kurtosis_average": moment(norm_abundance_average_list, moment=4),
                    "peak_average": max(norm_abundance_average, key=norm_abundance_average.get),
                    "max_fold_change_average": max_fold_change,
                    "residuals_average": residuals,
                    "R_squared_average": r_squared,
                    "residuals_all": residuals_all,
                    "R_squared_all": r_squared_all,
                    "curve_fold_change": curve_fold_change,
                    "curve_peak": curve_peak,
                }

        except Exception as e:
            print("+++++ EXCEPTION")
            print(e)
            # print(uniprot_accession)
            # print(norm_abundance_average)
            # print(norm_abundance_average_list)
            # print(range(0, len(norm_abundance_average)))

        return metrics
