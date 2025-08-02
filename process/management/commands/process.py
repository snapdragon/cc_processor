###########################################################################
# With thanks to ifigenia-t for the inspiration and a variety of functions.
# https://github.com/ifigenia-t
###########################################################################

import json
import logging
import math
import statistics
from time import sleep
import io
from collections import defaultdict

import numpy as np
import pandas as pd
import requests
from django.core.management.base import BaseCommand, CommandError
from django.db.models.query import QuerySet
from scipy import stats
from scipy.signal import periodogram
from scipy.special import comb
from scipy.stats import f_oneway, moment
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler

from process.constants import (
    ABUNDANCES_IMPUTED,
    ABUNDANCES_NORMALISED_LOG2_ARREST,
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ABUNDANCES_NORMALISED_MEDIAN,
    ABUNDANCES_NORMALISED_MIN_MAX,
    ABUNDANCES_NORMALISED_ZERO_MAX,
    ABUNDANCES_RAW,
    ANOVA,
    CURVE_FOLD_CHANGE,
    DEFAULT_FISHER_STATS,
    F_STATISTICS,
    FISHER_G,
    FREQUENCY,
    G_STATISTIC,
    ICR_ABUNDANCE_REP_1,
    ICR_ABUNDANCE_REP_2,
    P_VALUE,
    PHOSPHO_MAX_Q,
    PHOSPHO_MAX_Q_SL,
    PHOSPHO_MEDIAN,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN,
    PROTEIN_MAX_Q,
    PROTEIN_MAX_Q_SL,
    PROTEIN_MEDIAN,
    PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
    PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX,
    Q_VALUE,
    CCD_BATCH_SIZE,
    PROJECT_SL,
)
from process.models import (
    Abundance,
    Phospho,
    Project,
    Protein,
    Replicate,
    SampleStage,
    Statistic,
    StatisticType,
    UniprotData,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

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
            action="store_true",
        )
        parser.add_argument(
            "--calculate-proteins",
            help="Calculate and store protein results",
            action="store_true",
        )
        parser.add_argument(
            "--calculate-phosphos",
            help="Calculate and store phospho results",
            action="store_true",
        )
        parser.add_argument(
            "--calculate-phospho-medians",
            help="Calculate and store phospho medians",
            action="store_true",
        )
        parser.add_argument(
            "--calculate-batch",
            help="Calculate and store batch values",
            action="store_true",
        )
        parser.add_argument(
            "--run-all",
            help="Run all calculations and fetch all references",
            action="store_true",
        )
        parser.add_argument(
            "--fetch-references", help="Fetch references", action="store_true"
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
        run_all = options["run_all"]
        fetch_references = options["fetch_references"]

        if with_bugs and project_name != "ICR":
            raise CommandError("Only an ICR project can run --with-bugs")

        if not calculate_protein_medians and not calculate_proteins and not calculate_phospho_medians and not calculate_phosphos and not calculate_batch and not run_all and not fetch_references:
            raise CommandError("Add at least one flag to perform an action (typically --run-all)")

        project = Project.objects.get(name=project_name)

        if not project.processable:
            if (
                calculate_protein_medians
                or calculate_proteins
                or calculate_phospho_medians
                or calculate_phosphos
                or calculate_batch
                or run_all
            ):
                raise Exception(
                    f"Project {project.name} is non-processable, only --fetch-references can be run against it."
                )

        logger.info(f"Processing for project {project_name}, with bugs {with_bugs}")

        replicates = Replicate.objects.filter(project=project, mean=False)
        sample_stages = SampleStage.objects.filter(project=project).order_by("rank")

        project.with_bugs = with_bugs
        project.save()

        if calculate_protein_medians or run_all:
            self._calculate_medians(project, replicates, sample_stages, True, with_bugs)

        # N.B. there must be protein medians for this to run
        if calculate_proteins or run_all:
            if accession_number:
                logger.info(f"Processing protein {accession_number}")

                protein = Protein.objects.get(
                    project=project, accession_number=accession_number
                )

                self._calculate_abundances_metrics(
                    replicates,
                    sample_stages,
                    protein=protein,
                    phospho=None,
                    with_bugs=with_bugs,
                )
            else:
                logger.info("Processing all proteins")

                proteins = Protein.objects.filter(project=project).iterator(
                    chunk_size=100
                )

                for i, protein in enumerate(proteins):
                    if not i % 1000:
                        logger.info(
                            f"Calculating protein {i} {protein.accession_number}"
                        )

                    self._calculate_abundances_metrics(
                        replicates,
                        sample_stages,
                        protein=protein,
                        phospho=None,
                        with_bugs=with_bugs,
                    )

        if calculate_phospho_medians or run_all:
            self._calculate_medians(
                project, replicates, sample_stages, False, with_bugs
            )

        if calculate_phosphos or run_all:
            if accession_number:
                logger.info(f"Processing phospho protein {accession_number}")

                protein = Protein.objects.get(
                    project=project, accession_number=accession_number
                )

                self._calculate_phosphos(
                    replicates,
                    sample_stages,
                    protein,
                    with_bugs,
                )

            else:
                logger.info("Processing all phosphos")

                proteins = Protein.objects.filter(project=project).iterator(
                    chunk_size=100
                )

                for i, protein in enumerate(proteins):
                    if not i % 1000:
                        logger.info(
                            f"Calculating phospho protein {i} {protein.accession_number}"
                        )

                    self._calculate_phosphos(
                        replicates,
                        sample_stages,
                        protein,
                        with_bugs,
                    )

        if calculate_batch or run_all:
            self._calculate_protein_q_and_fisher_g(project, replicates, sample_stages)

            self._calculate_phospho_q_and_fisher_g(project, replicates, sample_stages)

            self._add_oscillations(project, replicates, sample_stages, with_bugs)

        if run_all or fetch_references:
            # Hits UniProt a lot, only run if necessary
            # self._get_uniprot_data(project)

            self._find_CCDs(project, True)
            self._find_CCDs(project, False)

            self._find_mitotic_cell_cycle(project, True)
            self._find_mitotic_cell_cycle(project, False)

            # self._fetch_corum_complexes(project)

            # Hits UniProt a lot, only run if necessary
            # # self._fetch_GO_locations(project)

            self._fetch_GO_enrichment(project, True)
            self._fetch_GO_enrichment(project, False)

            # self._generate_pcas(project)


    # # TODO - now taken care of by _get_uniprot_data
    # TODO - not really, _get_uniprot_data needs to be changed
    # def _fetch_GO_locations(self, project):
    #     Protein.objects.filter(project=project).update(GO_locations=None)

    #     proteins = Protein.objects.filter(
    #         project=project,
    #         ccd=True
    #     )

    #     num_proteins = 0

    #     for protein in proteins:
    #         if not num_proteins % 100:
    #             print(f"Fetching {num_proteins}")

    #         num_proteins += 1

    #         url = f"https://rest.uniprot.org/uniprotkb/{protein.accession_number}.json"
    #         response = requests.get(url)
    #         if response.status_code != 200:
    #             print(f"Error fetching {protein.accession_number}")
    #             return []
    #         data = response.json()
    
    #         go_terms = []

    #         for db_ref in data.get("uniProtKBCrossReferences", []):
    #             if db_ref["database"] == "GO":
    #                 go_id = db_ref["id"]

    #                 term_info = db_ref["properties"][0]["value"]  # e.g., "C:cytoplasm"
        
    #                 term_type, term_name = term_info.split(":", 1)

    #                 if term_type == "C":
    #                     go_terms.append((go_id, term_name))

    #         protein.GO_locations = go_terms
    #         protein.save()


    def _find_CCDs(self, project, is_protein):
        if is_protein:
            Protein.objects.filter(project=project).update(ccd=False)

            # SL has three replicates so a lower threshold is needed
            if project.name == PROJECT_SL:
                max_q = PROTEIN_MAX_Q_SL
            else:
                max_q = PROTEIN_MAX_Q

            statistics = Statistic.objects.filter(
                statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                protein__project=project,
            )
        else:
            Phospho.objects.filter(protein__project=project).update(ccd=False)

            if project.name == PROJECT_SL:
                max_q = PHOSPHO_MAX_Q_SL
            else:
                max_q = PHOSPHO_MAX_Q

            statistics = Statistic.objects.filter(
                statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                phospho__protein__project=project,
            )

        num_stats = 0

        print(f"MAX Q: {max_q}")

        for statistic in statistics:
            if not num_stats % 1000:
                logger.info(f"Generating relevance {num_stats}")

            num_stats += 1

            # Apparently q_values are rounded down, so no 'round' here
            if (
                statistic.metrics is None
                or statistic.metrics.get(ANOVA) is None
                or statistic.metrics[ANOVA].get(Q_VALUE) is None
                or statistic.metrics[ANOVA][Q_VALUE] >= max_q
            ):
                continue

            if (
                statistic.metrics.get(CURVE_FOLD_CHANGE) is None
                or self._round(statistic.metrics[CURVE_FOLD_CHANGE], 1) < 2
            ):
                continue

            if is_protein:
                statistic.protein.ccd = True
                statistic.protein.save()
            else:
                statistic.phospho.ccd = True
                statistic.phospho.save()

    def _find_mitotic_cell_cycle(self, project, is_protein):
        if is_protein:
            logger.info('Finding proteins with mitotic cell cycle GO term')
        else:
            logger.info('Finding phospho proteins with mitotic cell cycle GO term')

        matched_accessions = []

        if is_protein:
            Protein.objects.filter(project=project).update(mitotic_cell_cycle=False)

            records = Protein.objects.filter(
                project = project,
                ccd = True,
            )
        else:
            Protein.objects.filter(project=project).update(phospho_mitotic_cell_cycle=False)

            records = Phospho.objects.filter(
                protein__project = project,
                ccd = True,
            )

        for i in range(0, len(records), CCD_BATCH_SIZE):
            batch = records[i:i + CCD_BATCH_SIZE]

            if is_protein:
                accession_numbers = {p.accession_number for p in batch}
            else:
                accession_numbers = {p.protein.accession_number for p in batch}

            cleaned_accession_numbers = []

            # TODO - this should be used in the other place there is a check for
            #   a ;
            for an in list(accession_numbers):
                if ';' in an:
                    # Get first accession number only
                    cleaned_accession_numbers.extend(an.split(";"))
                else:
                    cleaned_accession_numbers.append(an)

            query = " OR ".join(f"accession:{an}" for an in cleaned_accession_numbers)

            url = "https://rest.uniprot.org/uniprotkb/search"

            params = {
                "query": query,
                "fields": "accession,go",
                "format": "tsv",
                "size": len(cleaned_accession_numbers)
            }

            response = requests.get(url, params=params)
            response.raise_for_status()
            df = pd.read_csv(io.StringIO(response.text), sep="\t")

            # Filter rows where GO term description contains 'mitotic cell cycle'
            if 'Gene Ontology (GO)' in df.columns:
                for _, row in df.iterrows():
                    go_term = str(row['Gene Ontology (GO)']).lower()

                    if 'mitotic cell cycle' in go_term:
                        accession_number = row['Entry']

                        try:
                            protein = Protein.objects.get(
                                project=project,
                                accession_number__contains = accession_number,
                            )

                            if is_protein:
                                protein.mitotic_cell_cycle = True
                            else:
                                protein.phospho_mitotic_cell_cycle = True

                            protein.save()
                        except Exception:
                            print(f"Can't get protein for accession_number {accession_number}")
            else:
                print("+++ NOT FOUND")



    def _fetch_corum_complexes(self, project):
        Protein.objects.filter(project=project).update(corum_results=None)

        # TODO - need phospho protein lists too?
        proteins = Protein.objects.filter(
            project=project,
            
        )

        accession_numbers = [p.accession_number for p in proteins]

        try:
            df = pd.read_csv(
                "data/corum_humanComplexes.txt", sep="\t", low_memory=False
            )
            df["subunits_uniprot_id"] = df["subunits_uniprot_id"].fillna("")
        except FileNotFoundError:
            print("CORUM file not found.")
            return None

        corum_results = {accession_number: [] for accession_number in accession_numbers}
        for _, row in df.iterrows():
            members = row["subunits_uniprot_id"].split(";")
            for accession_number in accession_numbers:
                if accession_number in members:
                    corum_results[accession_number].append((row["complex_id"], row["complex_name"]))

        for accession_number, complexes in corum_results.items():
            if not complexes:
                continue

            protein = Protein.objects.get(
                project=project,
                accession_number = accession_number
            )

            protein.corum_results = complexes
            protein.save()





    def _fetch_GO_enrichment(self, project, is_protein):
        if is_protein:
            logger.info(f"Fetching protein GO enrichment for project {project.name}")
        else:
            logger.info(
                f"Fetching phospho protein GO enrichment for project {project.name}"
            )

        if is_protein:
            proteins = Protein.objects.filter(
                project=project,
                ccd = True
            )

            accession_numbers = [p.accession_number for p in proteins]
        else:
            phosphos = Phospho.objects.filter(
                protein__project = project,
                ccd = True
            )

            accession_numbers = [p.protein.accession_number for p in phosphos]

        if not len(accession_numbers):
            logger.info(
                f"Project {project.name} has an empty accession number list for GO enrichment"
            )
            return

        genes = []

        for accession_number in accession_numbers:
            try:
                gene = self._get_gene_name(accession_number)

                genes.append(gene)
            except Exception:
                # Some proteins just don't have a gene in uniprot, skip
                logger.info(f"Can't find uniprot for {accession_number}")
                continue

        if not len(genes):
            logger.info(
                f"Project {project.name} can't find any genes for GO enrichment"
            )
            return

        payload = {
            "organism": "hsapiens",
            "query": genes,
            "sources": ["GO:BP", "GO:MF", "GO:CC", "Reactome"],
        }

        try:
            response = requests.post(
                "https://biit.cs.ut.ee/gprofiler/api/gost/profile/", json=payload
            )
            response.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"g:Profiler API request failed: {e}")
            return

        data = response.json()

        output = []

        for result in data.get("result", []):
            source = result.get("source")
            name = result.get("name")
            p_value = result.get("p_value")

            output.append(
                {
                    "source": source,
                    "name": name,
                    "p_value": p_value,
                }
            )

        if is_protein:
            project.protein_go_list = output
        else:
            project.phospho_protein_go_list = output

        project.save()

    def _generate_pcas(self, project):
        abundances = Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
            statistic__protein__project = project,
            replicate__mean = False,
            reading__isnull=False,
            # statistic__protein__mitotic_cell_cycle = True,
            # statistic__protein__ccd = True,
            # statistic__protein__accession_number__istartswith="P0"
        ).order_by("replicate__id", "sample_stage__rank").iterator(chunk_size=100)

        self._pca(project, abundances, True)

        abundances = Abundance.objects.filter(
            statistic__statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
            statistic__phospho__protein__project = project,
            replicate__mean = False,
            reading__isnull=False,
            # statistic__protein__mitotic_cell_cycle = True,
            # statistic__protein__ccd = True,
            # statistic__protein__accession_number__istartswith="P0"
        ).order_by("replicate__id", "sample_stage__rank")

        self._pca(project, abundances, False)


    def _pca(self, project, abundances, is_protein):
        abundance_table = {}

        for abundance in abundances:
            if is_protein:
                pan = abundance.statistic.protein.accession_number
            else:
                pan = abundance.statistic.phospho.mod

            reading = abundance.reading

            if abundance_table.get(pan) is None:
                abundance_table[pan] = {}

            abundance_table[pan][f"{abundance.replicate.name} {abundance.sample_stage.name}"] = reading

        df = pd.DataFrame(abundance_table)
        df = df.T


        if False:
            # Code for imputing values
            print("Dimensions before imputation:", df.shape)
            print("Missing values per column:", df.isnull().sum().sum())

            # # Impute missing values using median strategy
            # # You can change strategy to 'mean', 'most_frequent', or 'constant'
            # imputer = SimpleImputer(strategy='most_frequent')
            # df_imputed = pd.DataFrame(
            #     imputer.fit_transform(df),
            #     index=df.index,
            #     columns=df.columns
            # )

            from sklearn.impute import KNNImputer
    
            imputer = KNNImputer(n_neighbors=5)
            df_imputed = pd.DataFrame(
                imputer.fit_transform(df),
                index=df.index,
                columns=df.columns
            )

            print("Dimensions after imputation:", df_imputed.shape)
            print("Missing values after imputation:", df_imputed.isnull().sum().sum())

            df_T = df_imputed.T
        else:
            df = df.dropna()

            print("Dimensions of PCA dataframe:", df.shape)

            df_T = df.T


        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(df_T)

        sample_names = df_T.index.tolist()
        samples_info = []

        if project.name == "SL":
            num_stages = 10
        else:
            num_stages = 8

        for name in sample_names:
            idx = sample_names.index(name)

            if project.name == "SL":
                short_name = name.removeprefix("One ")
                short_name = short_name.removeprefix("Two ")
                short_name = short_name.removeprefix("Three ")
            else:
                short_name = name.removeprefix("abundance_rep_1")
                short_name = short_name.removeprefix("abundance_rep_2")

            if idx < num_stages:
                rep = 1
            elif idx < (num_stages * 2):
                rep = 2
            else:
                rep = 3

            stage = (idx - ( (rep-1) * num_stages)) + 1

            samples_info.append({
                "Sample": name,
                "PC1": pca_result[idx, 0],
                "PC2": pca_result[idx, 1],
                "Rep": f"{rep}",
                "Stage": f"{stage}",
                "Label": f"{short_name}"  # Optional: convert to time labels like 0h, 2h, etc.
            })

        context = {
            "pca_data_json": samples_info,
            "pc1_var": f"{pca.explained_variance_ratio_[0]*100:.1f}%",
            "pc2_var": f"{pca.explained_variance_ratio_[1]*100:.1f}%"
        }

        if is_protein:
            project.protein_pca = context
        else:
            project.phospho_pca = context
        
        project.save()

    def _calculate_phosphos(self, replicates, sample_stages, protein, with_bugs):
        phosphos = Phospho.objects.filter(protein=protein)

        for phospho in phosphos:
            self._calculate_abundances_metrics(
                replicates, sample_stages, None, phospho, with_bugs
            )

    # TODO - similar to other functionality, consolidate
    def _calculate_phospho_q_and_fisher_g(self, project, replicates, sample_stages):
        logger.info("Calculating phospho q and fisher G values")
        # TODO - blank all q_values in DB?

        fisher_g_stats = self._calculate_fisher_g(
            project, replicates, sample_stages, phospho=True
        )

        # Get all phospho abundances for this project for log2 mean
        statistics = Statistic.objects.filter(
            statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
            phospho__protein__project=project,
        )

        anova_stats: dict = {}

        for statistic in statistics:
            site_key = self._get_site_key(statistic)

            anova_stats[site_key] = {P_VALUE: statistic.metrics[ANOVA][P_VALUE]}

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
            statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein__project=project,
        )

        fisher_g_stats = self._calculate_fisher_g(project, replicates, sample_stages)

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

    # TODO - this can probably be refactored. All it does is get the readings in one
    #   form and then change them to another. Might as well go straight to the other.
    def _create_abundance_dataframe(
        self, project, replicates, sample_stages, phospho, phospho_ab, phospho_reg
    ):
        abundance_table = {}

        # TODO - could combine these two loops using a function to generate the keys
        if phospho:
            # Determine which phospho readings to use
            statistic_type_name = ABUNDANCES_NORMALISED_LOG2_MEAN

            if phospho_ab:
                statistic_type_name = PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN
            elif phospho_reg:
                statistic_type_name = (
                    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN
                )

            # Get all log2 mean abundances for all chosen phospho readings for this project
            abundances = Abundance.objects.filter(
                statistic__statistic_type__name=statistic_type_name,
                statistic__phospho__protein__project=project,
                replicate__mean = False,
                reading__isnull=False,
            ).iterator(chunk_size=100)

            num_abundances = 0

            for abundance in abundances:
                if not num_abundances % 100000:
                    logger.info(
                        f"Processing dataframe phospho abundance {num_abundances}"
                    )

                num_abundances += 1

                site_key = self._get_site_key(abundance.statistic)

                if abundance_table.get(site_key) is None:
                    abundance_table[site_key] = {}

                rep_stage_name = self._rep_stage_name(abundance)

                abundance_table[site_key][rep_stage_name] = abundance.reading
        else:
            # Get all log2 mean abundances for all proteins for this project
            abundances = Abundance.objects.filter(
                statistic__statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
                statistic__protein__project=project,
                replicate__mean = False,
                reading__isnull = False
            ).iterator(chunk_size=100)

            num_abundances = 0

            for abundance in abundances:
                if not num_abundances % 100000:
                    logger.info(
                        f"Processing dataframe protein abundance {num_abundances}"
                    )

                num_abundances += 1

                pan = abundance.statistic.protein.accession_number

                # # Not all proteins have protein results, some are phospho only
                # # TODO - is this necessary? Would phospho-only results have a record?
                # # TODO - can remove this by putting value__isnull=False in the query
                # if abundance.reading is None:
                #     continue

                if abundance_table.get(pan) is None:
                    abundance_table[pan] = {}

                rep_stage_name = self._rep_stage_name(abundance)

                abundance_table[pan][rep_stage_name] = abundance.reading

        abundance_table_df = pd.DataFrame(abundance_table)
        abundance_table_df = abundance_table_df.T

        new_cols = []

        # TODO - is this ordering correct?
        for sample_stage in sample_stages:
            for replicate in replicates:
                new_cols.append(f"{replicate.name}_{sample_stage.name}")

        # Rearrange the column order so replicates are near eatch other
        try:
            abundance_table_df = abundance_table_df[new_cols]
        except Exception as e:
            logger.error("DF FAILED")
            logger.error(abundance_table_df)
            logger.error(e)
            exit()

        return abundance_table_df

    # # TODO - is this really needed any more?
    # def _tp(self, abundances):
    #     """
    #     Creates a list of each abundance of a stage name across replicates
    #     """
    #     res = []

    #     for abundance in abundances:
    #         res.append(abundance.reading)

    #     return res

    def _calculate_ANOVA(
        self, statistic_type_name, sample_stages, protein=None, phospho=None
    ):
        # Defaults if not enough replicates
        p_value = 1
        f_statistic = 1

        stat_log2_mean = self._get_statistic(
            statistic_type_name,
            protein=protein,
            phospho=phospho,
        )

        abundances = Abundance.objects.filter(
            statistic=stat_log2_mean,
            replicate__mean=False,
            reading__isnull=False,
        ).order_by("sample_stage__rank")

        # TODO - not a good name
        stage_names = []

        try:
            for sample_stage in sample_stages:
                abundances.filter(sample_stage=sample_stage),

                stage_names.append([a.reading for a in abundances.filter(sample_stage=sample_stage)])

                # stage_names.append(
                #     self._tp(
                #         abundances.filter(sample_stage=sample_stage),
                #     )
                # )

            # Each entry must have at least two points for f_oneway to work
            timepoints = [x for x in stage_names if x != [] and len(x) > 1]

            if len(timepoints) > 1:
                one_way_anova = f_oneway(*timepoints)

                f_statistic = one_way_anova[0].item()
                p_value = one_way_anova[1].item()

                if np.isnan(p_value):
                    p_value = 1

            metrics = stat_log2_mean.metrics
            metrics[ANOVA] = {P_VALUE: p_value, F_STATISTICS: f_statistic}

            stat_log2_mean.save()

        except Exception as e:
            logger.error("Error in _calculate_ANOVA")
            logger.error(e)

    def _calculate_metrics(
        self,
        statistic_type_name: str,
        replicates,
        sample_stages,
        protein=None,
        phospho=None,
    ):
        metrics = {}

        if protein:
            stat = self._get_statistic(statistic_type_name, protein=protein)
        else:
            stat = self._get_statistic(statistic_type_name, phospho=phospho)

        abundances_all = Abundance.objects.filter(statistic=stat).order_by(
            "replicate__id", "sample_stage__rank"
        )

        # Clear previous metrics
        stat.metrics = {}
        stat.save()

        abundances_non_mean = abundances_all.filter(replicate__mean=False)

        abundance_averages = {
            a.sample_stage.name: a.reading
            for a in abundances_all.filter(replicate__mean=True)
            if a.reading is not None
        }
        abundance_averages_list = [val for val in abundance_averages.values()]

        if len(abundances_non_mean.values("replicate__id").distinct()) < 2:
            # There aren't as least two replicates to calculate, give up
            return

        if not len(abundance_averages_list) or not len(abundances_non_mean):
            # Nothing to work with
            return

        std = statistics.stdev(
            [a.reading for a in abundances_non_mean if a.reading is not None]
        )

        try:
            polyfit_result = np.polyfit(
                range(0, len(abundance_averages)),
                abundance_averages_list,
                2,
                full=True,
            )

            if len(polyfit_result[1]):
                residuals = polyfit_result[1][0]
            else:
                # TODO - is this a meaningful value? It's just the mean.
                residuals = int(len(sample_stages) / 2)

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

            curve_fold_change, curve_peak = self._calculate_curve_fold_change(
                abundances_non_mean,
                replicates,
            )

            residuals_all, r_squared_all = self._calculate_residuals_R2_all(
                abundances_non_mean,
                replicates,
            )

            metrics = {
                "standard_deviation": std,
                "variance_average": round(moment(abundance_averages_list, moment=2), 2),
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
    def _calculate_residuals_R2_all(self, abundances, replicates):
        residuals_all = None
        r_squared_all = None

        x, y, _ = self._generate_xs_ys(abundances, replicates)

        if len(x) != len(y):
            return residuals_all, r_squared_all

        p = np.poly1d(np.polyfit(x, y, 2))
        curve_abundances = p(x)
        residuals_all = np.polyfit(x, y, 2, full=True)[1][0]
        r_squared_all = round(r2_score(y, curve_abundances), 2)

        return residuals_all.item(), r_squared_all

    # TODO - this is straight up lifted from ICR, change and ideally use a library
    def _calculate_curve_fold_change(self, abundances, replicates):
        curve_fold_change = None
        curve_peak = None

        x, y, stage_names_map = self._generate_xs_ys(
            abundances,
            replicates,
        )

        if len(x) == len(y):
            p = np.poly1d(np.polyfit(x, y, 2))
            curve_abundances = p(x)

            # find the timepoint peak of the curve
            curve_index = x[list(curve_abundances).index(max(curve_abundances))]
            for time_point, index in stage_names_map.items():
                if index == curve_index:
                    curve_peak = time_point

            # Calculate the fold change from the curve
            curve_fold_change = max(curve_abundances) / max(0.05, min(curve_abundances))
            curve_fold_change = curve_fold_change.item()

        return curve_fold_change, curve_peak

    # TODO - this is straight up lifted from ICR. Replace it, ideally with a library call
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
        replicates: QuerySet[Replicate],
        protein=None,
        phospho=None,
    ):
        _, stat_imputed = self._clear_and_fetch_stats(
            ABUNDANCES_IMPUTED, protein=protein, phospho=phospho
        )

        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MIN_MAX,
            protein=protein,
            phospho=phospho,
        ).order_by("sample_stage__rank")

        abundances_by_rep: dict = {}

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
                        #   Make this a 'with_bugs'?
                        next_abundance = abs[j % len(abs)]

                        if next_abundance.reading is not None:
                            next = (j, next_abundance.reading)
                            break

                    if last and next:
                        # Linear imputation between nearest timepoints
                        last_offset, last_reading = last
                        next_offset, next_reading = next
                        step_height = (last_reading - next_reading) / (
                            last_offset + next_offset
                        )
                        reading = next_offset * step_height + next_reading

                Abundance.objects.create(
                    statistic=stat_imputed,
                    replicate=abundance.replicate,
                    sample_stage=abundance.sample_stage,
                    reading=reading,
                )

    def _calculate_zero_or_min_normalisation(
        self, replicates, protein=None, phospho=None, zero_min=False
    ):
        if zero_min:
            stat_type_cf = ABUNDANCES_NORMALISED_ZERO_MAX
            stat_type_get = ABUNDANCES_NORMALISED_MEDIAN
        else:
            stat_type_cf = ABUNDANCES_NORMALISED_MIN_MAX
            stat_type_get = ABUNDANCES_NORMALISED_LOG2_MEAN

        _, stat = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_ZERO_MAX, protein=protein, phospho=phospho
        )

        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MEDIAN, protein=protein, phospho=phospho
        )

        # # TODO - can be tidied, only constants are different
        # if zero_min:
        #     _, stat = self._clear_and_fetch_stats(
        #         ABUNDANCES_NORMALISED_ZERO_MAX, protein=protein, phospho=phospho
        #     )

        #     abundances = self._get_abundances(
        #         ABUNDANCES_NORMALISED_MEDIAN, protein=protein, phospho=phospho
        #     )
        # else:
        #     _, stat = self._clear_and_fetch_stats(
        #         ABUNDANCES_NORMALISED_MIN_MAX, protein=protein, phospho=phospho
        #     )

        #     abundances = self._get_abundances(
        #         ABUNDANCES_NORMALISED_LOG2_MEAN, protein=protein, phospho=phospho
        #     )

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
                    reading_normalised = 0.5
                else:
                    # TODO - why rounded?
                    reading_normalised = round((reading - min_value) / denominator, 4)

                Abundance.objects.create(
                    statistic=stat,
                    replicate=ab.replicate,
                    sample_stage=ab.sample_stage,
                    reading=reading_normalised,
                )

    def _calculate_relative_log2_normalisation(self, protein=None, phospho=None):
        _, stat_normalised_log2_mean = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_LOG2_MEAN, protein=protein, phospho=phospho
        )

        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MEDIAN,
            protein=protein,
            phospho=phospho,
        )

        total_abundances = 0
        total_lengths = 0

        for abundance in abundances:
            if abundance.reading is not None and abundance.reading != 0:
                total_abundances += math.log2(abundance.reading)
                total_lengths += 1

        mean = None

        if total_lengths != 0:
            mean = total_abundances / total_lengths

        for abundance in abundances:
            if abundance.reading is None or abundance.reading == 0:
                continue

            normalised_abundance = self._round(math.log2(abundance.reading) - mean)

            Abundance.objects.create(
                statistic=stat_normalised_log2_mean,
                replicate=abundance.replicate,
                sample_stage=abundance.sample_stage,
                reading=normalised_abundance,
            )

    def _calculate_arrest_log2_normalisation(self, protein=None, phospho=None):
        ARRESTING_AGENT = "Nocodozole"

        _, stat_normalised_log2_arrest = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_LOG2_ARREST, protein=protein, phospho=phospho
        )

        if protein:
            # TODO - this is a hack, add the field to the Project model.
            if protein.project.name == "ICR":
                ARRESTING_AGENT = "Palbo"
        else:
            if phospho.protein.project.name == "ICR":
                ARRESTING_AGENT = "Palbo"

        abundances = self._get_abundances(
            ABUNDANCES_NORMALISED_MEDIAN, protein=protein, phospho=phospho
        )

        for abundance in abundances:
            reading = abundance.reading

            # Get the arresting agent abundance for this replicate
            arrest_abundance = abundances.filter(
                sample_stage__name=ARRESTING_AGENT, replicate=abundance.replicate
            ).first()

            # Not all replicates have an arresting agent value
            # TODO - do we need to check for the other two Nones?
            if (
                reading is not None
                and arrest_abundance is not None
                and arrest_abundance.reading is not None
                and reading != 0
                and arrest_abundance.reading != 0
            ):
                log2_reading = self._round(
                    math.log2(reading / arrest_abundance.reading)
                )

                Abundance.objects.create(
                    statistic=stat_normalised_log2_arrest,
                    replicate=abundance.replicate,
                    sample_stage=abundance.sample_stage,
                    reading=log2_reading,
                )

    def _get_statistic_type(self, name):
        return StatisticType.objects.get(name=name)

    def _get_statistic(
        self, statistic_type_name, project=None, protein=None, phospho=None
    ):
        statistic_type = self._get_statistic_type(statistic_type_name)

        if project:
            stat, _ = Statistic.objects.get_or_create(
                statistic_type=statistic_type, project=project
            )
        elif protein:
            stat, _ = Statistic.objects.get_or_create(
                statistic_type=statistic_type, protein=protein
            )
        elif phospho:
            stat, _ = Statistic.objects.get_or_create(
                statistic_type=statistic_type, phospho=phospho
            )
        else:
            raise Exception("_get_statistic needs a project, protein or phospho")

        return stat

    def _get_abundances(
        self, statistic_type_name, project=None, protein=None, phospho=None
    ):
        statistic = self._get_statistic(statistic_type_name, project, protein, phospho)

        return Abundance.objects.filter(
            statistic=statistic,
            reading__isnull = False
        )

    def _calculate_normalised_medians(self, protein=None, phospho=None):
        _, stat_normalised_medians = self._clear_and_fetch_stats(
            ABUNDANCES_NORMALISED_MEDIAN, protein=protein, phospho=phospho
        )

        if protein:
            stat_medians = self._get_statistic(PROTEIN_MEDIAN, project=protein.project)
        else:
            stat_medians = self._get_statistic(
                PHOSPHO_MEDIAN, project=phospho.protein.project
            )

        readings = self._get_abundances(
            ABUNDANCES_RAW, protein=protein, phospho=phospho
        )

        for prr in readings:
            reading = prr.reading

            # TODO - is this if statement needed?
            if reading is not None:
                median = Abundance.objects.get(
                    statistic=stat_medians,
                    replicate=prr.replicate,
                    sample_stage=prr.sample_stage,
                )

                normalised_reading = reading / median.reading

                Abundance.objects.create(
                    statistic=stat_normalised_medians,
                    replicate=prr.replicate,
                    sample_stage=prr.sample_stage,
                    reading=normalised_reading,
                )

    def _calculate_means(
        self,
        statistic_type_name: str,
        protein=None,
        phospho=None,
        with_bugs: bool = False,
        imputed: bool = False,
    ):
        readings: dict = {}

        abundances = self._get_abundances(
            statistic_type_name, protein=protein, phospho=phospho
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
            mean_replicate = Replicate.objects.get(
                project=phospho.protein.project, mean=True
            )

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
                reading=mean,
            )

    def _generate_xs_ys(self, abundances, replicates):
        # TODO - why do we use the second replicate?
        abundances = abundances.filter(replicate=replicates[1])

        stage_names_map = {}

        x = []
        y = []
        for i, abundance in enumerate(abundances):
            stage_names_map[abundance.sample_stage.name] = i

            x.append(stage_names_map[abundance.sample_stage.name])
            if abundance.reading is not None:
                y.append(abundance.reading)
        x.sort()

        return x, y, stage_names_map

    def _add_oscillations(self, project, replicates, sample_stages, with_bugs):
        # TODO - check all logging statements for similarity to ICR
        logger.info("Adding Protein Oscillation Normalised Abundances")

        self._generate_protein_oscillation_metrics(
            project, replicates, sample_stages, with_bugs
        )

        fisher_g_stats = self._calculate_fisher_g(
            project, replicates, sample_stages, phospho=True, phospho_ab=True
        )

        ps_and_qs = {}

        num_proteins = 0

        statistics = Statistic.objects.filter(
            statistic_type__name=PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
            phospho__protein__project=project,
        ).iterator(chunk_size=100)

        for statistic in statistics:
            if not num_proteins % 1000:
                print(
                    f"Calculating protein oscillation {num_proteins} {statistic.phospho.protein.accession_number}"
                )

            num_proteins += 1

            phospho_key = self._get_site_key(statistic)

            if phospho_key not in ps_and_qs:
                ps_and_qs[phospho_key] = {}

            if (statistic.metrics.get(ANOVA) is not None) and (
                statistic.metrics[ANOVA].get(P_VALUE) is not None
            ):
                ps_and_qs[phospho_key][P_VALUE] = statistic.metrics[ANOVA][P_VALUE]

        ps_and_qs = self._add_q_value(ps_and_qs)

        num_proteins = 0

        statistics = Statistic.objects.filter(
            statistic_type__name=PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
            phospho__protein__project=project,
        ).iterator(chunk_size=100)

        # TODO - this code is very similar to others, consolidate?
        for statistic in statistics:
            if not num_proteins % 1000:
                print(
                    f"Calculating phospho oscillation {num_proteins} {statistic.phospho.protein.accession_number}"
                )

            num_proteins += 1

            phospho_key = self._get_site_key(statistic)

            q_value = 1

            if phospho_key in ps_and_qs:
                q_value = ps_and_qs[phospho_key][Q_VALUE]

            statistic.metrics[ANOVA][Q_VALUE] = q_value

            fisher_output = DEFAULT_FISHER_STATS

            if phospho_key in fisher_g_stats:
                fisher_output = fisher_g_stats[phospho_key]

            statistic.metrics[FISHER_G] = fisher_output

            statistic.save()

    def _calculate_protein_oscillation(
        self,
        protein,
        phospho,
        replicates,
        sample_stages,
        zero_or_log2,
        statistic_type_name,
    ):
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
                except Exception:
                    # logger.info("No protein or phospho for calculating oscillation")
                    continue

                oscillation = phospho_stage.reading - protein_stage.reading

                Abundance.objects.create(
                    statistic=stat,
                    replicate=replicate,
                    sample_stage=sample_stage,
                    reading=oscillation,
                )

    def _calculate_abundances_metrics(
        self, replicates, sample_stages, protein=None, phospho=None, with_bugs=False
    ):
        if not protein and not phospho:
            logger.error("Either a protein or a phospho must be provided.")

        if protein:
            Abundance.objects.filter(
                statistic__statistic_type__name=ABUNDANCES_RAW,
                statistic__protein=protein,
                replicate__mean=True,
            ).delete()
        else:
            Abundance.objects.filter(
                statistic__statistic_type__name=ABUNDANCES_RAW,
                statistic__phospho=phospho,
                replicate__mean=True,
            ).delete()

        self._calculate_normalised_medians(protein, phospho)

        self._calculate_arrest_log2_normalisation(protein, phospho)

        self._calculate_relative_log2_normalisation(protein, phospho)

        self._calculate_zero_or_min_normalisation(replicates, protein, phospho, False)

        self._calculate_zero_or_min_normalisation(replicates, protein, phospho, True)

        self._impute(
            replicates,
            protein,
            phospho,
        )

        self._calculate_means(
            ABUNDANCES_RAW, protein=protein, phospho=phospho, with_bugs=with_bugs
        )

        self._calculate_means(
            ABUNDANCES_NORMALISED_MEDIAN,
            protein=protein,
            phospho=phospho,
            with_bugs=with_bugs,
        )

        self._calculate_means(
            ABUNDANCES_NORMALISED_MIN_MAX,
            protein=protein,
            phospho=phospho,
            with_bugs=with_bugs,
        )

        self._calculate_means(
            ABUNDANCES_NORMALISED_ZERO_MAX,
            protein=protein,
            phospho=phospho,
            with_bugs=with_bugs,
        )

        self._calculate_means(
            ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein=protein,
            phospho=phospho,
            with_bugs=with_bugs,
        )

        self._calculate_means(
            ABUNDANCES_NORMALISED_LOG2_ARREST,
            protein=protein,
            phospho=phospho,
            with_bugs=with_bugs,
        )

        self._calculate_means(
            ABUNDANCES_IMPUTED,
            protein=protein,
            phospho=phospho,
            with_bugs=with_bugs,
            imputed=True,
        )

        self._calculate_metrics(
            ABUNDANCES_NORMALISED_LOG2_MEAN,
            replicates,
            sample_stages,
            protein,
            phospho,
        )

        self._calculate_metrics(
            ABUNDANCES_NORMALISED_ZERO_MAX,
            replicates,
            sample_stages,
            protein,
            phospho,
        )

        self._calculate_ANOVA(
            ABUNDANCES_NORMALISED_LOG2_MEAN, sample_stages, protein, phospho
        )

    def _generate_protein_oscillation_metrics(
        self, project, replicates, sample_stages, with_bugs
    ):
        phosphos = Phospho.objects.filter(protein__project=project).iterator(
            chunk_size=100
        )

        num_phosphos = 0

        for phospho in phosphos:
            if not num_phosphos % 1000:
                logger.info(
                    f"Processing oscillation {num_phosphos} {phospho.protein.accession_number} {phospho.mod}"
                )

            num_phosphos += 1

            for zero_or_log2, statistic_type_name in [
                [
                    ABUNDANCES_NORMALISED_ZERO_MAX,
                    PROTEIN_OSCILLATION_ABUNDANCES_ZERO_MAX,
                ],
                [
                    ABUNDANCES_NORMALISED_LOG2_MEAN,
                    PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
                ],
            ]:
                self._calculate_protein_oscillation(
                    phospho.protein,
                    phospho,
                    replicates,
                    sample_stages,
                    zero_or_log2,
                    statistic_type_name,
                )

                self._calculate_means(
                    statistic_type_name,
                    protein=None,
                    phospho=phospho,
                    with_bugs=with_bugs,
                )

                self._calculate_metrics(
                    statistic_type_name, replicates, sample_stages, None, phospho
                )

                self._calculate_ANOVA(statistic_type_name, sample_stages, None, phospho)

    def _add_q_value(self, data):
        df = pd.DataFrame(data)
        df = df.T

        df[Q_VALUE] = stats.false_discovery_control(df[P_VALUE])

        return df.to_dict("index")

    def _calculate_medians(
        self, project, replicates, sample_stages, is_protein, with_bugs
    ):
        if is_protein:
            logger.info("Calculating protein medians")
        else:
            logger.info("Calculating phospho medians")

        if is_protein:
            _, stat_prot_med = self._clear_and_fetch_stats(
                PROTEIN_MEDIAN, project=project
            )

            abundances = Abundance.objects.filter(
                statistic__protein__project=project,
                statistic__statistic_type__name=ABUNDANCES_RAW,
            ).iterator(chunk_size=100)
        else:
            _, stat_prot_med = self._clear_and_fetch_stats(
                PHOSPHO_MEDIAN, project=project
            )

            abundances = Abundance.objects.filter(
                statistic__phospho__protein__project=project,
                statistic__statistic_type__name=ABUNDANCES_RAW,
            ).iterator(chunk_size=100)

        rep_stage_abundances = {}

        for i, abundance in enumerate(abundances):
            if not i % 10000:
                if is_protein:
                    logger.info(f"Processing protein abundance for median {i}")
                else:
                    logger.info(f"Processing phospho abundance for median {i}")

            if rep_stage_abundances.get(abundance.replicate) is None:
                rep_stage_abundances[abundance.replicate] = {}

            if (
                rep_stage_abundances[abundance.replicate].get(abundance.sample_stage)
                is None
            ):
                rep_stage_abundances[abundance.replicate][abundance.sample_stage] = []

            rep_stage_abundances[abundance.replicate][abundance.sample_stage].append(
                abundance.reading
            )

        if is_protein and with_bugs:
            replicate1 = Replicate.objects.get(
                project=project, name=ICR_ABUNDANCE_REP_1
            )
            replicate2 = Replicate.objects.get(
                project=project, name=ICR_ABUNDANCE_REP_2
            )

            rep_stage_abundances[replicate1] = rep_stage_abundances[replicate2]

        for replicate in replicates:
            if rep_stage_abundances.get(replicate) is None:
                continue

            for sample_stage in sample_stages:
                if rep_stage_abundances[replicate].get(sample_stage) is None:
                    continue

                if not len(rep_stage_abundances[replicate][sample_stage]):
                    logger.error(
                        f"Median with no values (??) {replicate.name} {sample_stage.name}"
                    )
                    continue

                median = statistics.median(
                    rep_stage_abundances[replicate][sample_stage]
                )

                Abundance.objects.create(
                    statistic=stat_prot_med,
                    replicate=replicate,
                    sample_stage=sample_stage,
                    reading=median,
                )

    def _calculate_fisher_g(
        self,
        project,
        replicates,
        sample_stages,
        phospho=False,
        phospho_ab=False,
        phospho_reg=False,
    ):
        logger.info("Calculate Fisher G Statistics")

        abundance_fisher = self._create_abundance_dataframe(
            project, replicates, sample_stages, phospho, phospho_ab, phospho_reg
        )
        abundance_fisher = abundance_fisher.dropna()

        for index, row in abundance_fisher.iterrows():
            row_z = [i for i in row.tolist()]

            result = self._ptestg(row_z)

            abundance_fisher.loc[index, G_STATISTIC] = result["obsStat"]
            abundance_fisher.loc[index, P_VALUE] = result["pvalue"]
            abundance_fisher.loc[index, FREQUENCY] = result["freq"]

        q_value = stats.false_discovery_control(abundance_fisher[P_VALUE])

        abundance_fisher[Q_VALUE] = q_value

        cols = abundance_fisher.columns
        fisher_cols = [G_STATISTIC, P_VALUE, FREQUENCY, Q_VALUE]
        ab_col = [x for x in cols if x not in fisher_cols]
        abundance_fisher = abundance_fisher.drop(columns=ab_col)
        abundance_fisher_dict = abundance_fisher.to_dict("index")

        return abundance_fisher_dict

    def _get_uniprot_data(self, project):
        logger.info("Fetching uniprot data")

        proteins = Protein.objects.filter(project=project)

        for protein in proteins:
            try:
                # A problem with this is that, being accululative, it
                #   doesn't blank all Protein GO_locations values in advance.
                #   If it's a concern the table can be emptied.
                UniprotData.objects.get(accession_number=protein.accession_number)
            except Exception:
                logger.info(f"Fetching uniprot for {protein.accession_number}")

                url = f"https://rest.uniprot.org/uniprotkb/{protein.accession_number}.json"

                response = requests.get(url)

                if response.status_code == 200:
                    data = response.json()

                    protein_name = (
                        data.get("proteinDescription", {})
                        .get("recommendedName", {})
                        .get("fullName", {})
                        .get("value", "Unknown")
                    )
                    gene_name = (
                        data.get("genes", [{}])[0]
                        .get("geneName", {})
                        .get("value", "Unknown")
                        if data.get("genes")
                        else "Unknown"
                    )

                    UniprotData.objects.create(
                        accession_number=protein.accession_number,
                        gene_name=gene_name,
                        protein_name=protein_name,
                    )

                    go_terms = []

                    for db_ref in data.get("uniProtKBCrossReferences", []):
                        if db_ref["database"] == "GO":
                            go_id = db_ref["id"]

                            term_info = db_ref["properties"][0]["value"]
        
                            term_type, term_name = term_info.split(":", 1)

                            if term_type == "C":
                                go_terms.append((go_id, term_name))

                    protein.GO_locations = go_terms
                    protein.save()
                else:
                    logger.error(
                        f"Failed to retrieve data from UniProt for {protein.accession_number}. Status code: {response.status_code}"
                    )

                sleep(0.1)

    def _get_gene_name(self, accession_number):
        uniprot_data = UniprotData.objects.get(accession_number=accession_number)

        return uniprot_data.gene_name

    def _clear_and_fetch_stats(
        self, statistic_type_name, project=None, protein=None, phospho=None
    ):
        statistic_type = self._get_statistic_type(statistic_type_name)

        if project:
            stat, _ = Statistic.objects.get_or_create(
                statistic_type=statistic_type, project=project
            )
        elif protein:
            stat, _ = Statistic.objects.get_or_create(
                statistic_type=statistic_type, protein=protein
            )
        elif phospho:
            stat, _ = Statistic.objects.get_or_create(
                statistic_type=statistic_type, phospho=phospho
            )
        else:
            raise Exception(
                "_clear_and_fetch_stats needs a project, protein or phospho"
            )

        Abundance.objects.filter(statistic=stat).delete()

        return statistic_type, stat

    def _round(self, value, places=4):
        return round(value, places)

    def _get_site_key(self, statistic: Statistic):
        return f"{statistic.phospho.protein.accession_number}_{statistic.phospho.mod}"

    def _rep_stage_name(self, abundance):
        return f"{abundance.replicate.name}_{abundance.sample_stage.name}"

    def _fisher_g_test(self, z):
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
        freqs, pgram_vals = periodogram(z, scaling="spectrum")

        # Remove Nyquist frequency for even-length input
        if n % 2 == 0:
            pgram_vals = pgram_vals[:-1]

        max_index = np.argmax(pgram_vals)
        g_stat = pgram_vals[max_index] / np.sum(pgram_vals)

        # Compute p-value
        p = int(np.floor(1 / g_stat))
        i_vals = np.arange(1, p + 1)
        terms = (
            comb(m, i_vals) * (-1) ** (i_vals - 1) * (1 - i_vals * g_stat) ** (m - 1)
        )
        p_value = np.sum(terms)

        return np.array([g_stat, p_value, freqs[max_index]])

    def _ptestg(self, z):
        """
        Perform periodicity test on input data `z` using Fisher's g-test.
        """
        z = np.atleast_2d(z).T  # Ensure z is a column vector
        n, m = z.shape

        obs_stat = np.zeros(m)
        p_values = np.zeros(m)
        freqs = np.zeros(m)

        for i in range(m):
            result = self._fisher_g_test(z[:, i])
            obs_stat[i] = result[0]
            p_values[i] = result[1]
            freqs[i] = result[2]

        return {
            "obsStat": obs_stat,
            "pvalue": p_values,
            "freq": freqs,
            "class": "Htest",
        }
