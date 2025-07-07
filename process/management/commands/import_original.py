import logging

import pandas as pd
from django.core.management.base import BaseCommand
import ijson
from decimal import Decimal

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
    PROTEIN_ABUNDANCES_RAW,
    PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
    NORMALISED,
    MEDIAN,
    PROTEIN_ABUNDANCES_NORMALISED_LOG2_ARREST,
    PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
    LOG2_MEAN,
    PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX,
    MIN_MAX,
    PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX,
    ZERO_MAX,
    PROTEIN_ABUNDANCES_IMPUTED,
    IMPUTED,
    METRICS,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

PROJECT_NAME = "Original"

class Command(BaseCommand):
    help = "import the output from the original for comparison"

    def handle(self, *args, **options):
        logger.info(f"Importing original output")

        project = Project.objects.get(name=PROJECT_NAME)
        stats_type_rp = StatisticType.objects.get(name=PROTEIN_ABUNDANCES_RAW)

        file_path = f"data/ICR/TimeCourse_Full_info_full_indented.json"

        Abundance.objects.filter(statistic__project=project).delete()
        Abundance.objects.filter(statistic__protein__project=project).delete()
        Abundance.objects.filter(statistic__phospho__protein__project=project).delete()
        Statistic.objects.filter(project=project).delete()
        Statistic.objects.filter(protein__project=project).delete()
        Statistic.objects.filter(phospho__protein__project=project).delete()
        Protein.objects.filter(project=project).delete()
        Phospho.objects.filter(protein__project=project).delete()

        replicates = Replicate.objects.filter(project=project)
        sample_stages = SampleStage.objects.filter(project=project).order_by('rank')

        replicates_by_name = {}
        sample_stages_by_name = {}

        for replicate in replicates:
            replicates_by_name[replicate.name] = replicate

        for sample_stage in sample_stages:
            sample_stages_by_name[sample_stage.name] = sample_stage

        num_proteins = 0

        with open(file_path, 'r') as f:
            for gene_name, gene_data in ijson.kvitems(f, ''):
                num_proteins += 1

                if not num_proteins % 1000:
                    print(f"Importing original protein result {num_proteins} {gene_name}")

                protein = Protein.objects.create(
                    project=project, accession_number=gene_name, is_contaminant=False
                )

                if not gene_data.get(PROTEIN_ABUNDANCES):
                    print(f"No protein abundances for {gene_name}")
                    continue

                pa = gene_data[PROTEIN_ABUNDANCES]

                if not gene_data.get(METRICS):
                    print(f"No metrics for {gene_name}")
                    continue

                pm = gene_data[METRICS]

                self._import_protein_readings(
                    replicates_by_name,
                    sample_stages_by_name,
                    protein,
                    PROTEIN_ABUNDANCES_RAW,
                    pa[RAW],
                )

                # Apparently not all originals have normalised medians
                #   Phospho-only imports maybe?
                if pa.get(NORMALISED) and pa[NORMALISED].get(MEDIAN):
                    self._import_protein_readings(
                        replicates_by_name,
                        sample_stages_by_name,
                        protein,
                        PROTEIN_ABUNDANCES_NORMALISED_MEDIAN,
                        pa[NORMALISED][MEDIAN],
                    )
                else:
                    print(f"No normalised medians for protein {protein.accession_number}")

                if pa.get(NORMALISED) and pa[NORMALISED].get("log2_palbo"):
                    self._import_protein_readings(
                        replicates_by_name,
                        sample_stages_by_name,
                        protein,
                        PROTEIN_ABUNDANCES_NORMALISED_LOG2_ARREST,
                        pa[NORMALISED]["log2_palbo"],
                    )
                else:
                    print(f"No normalised log2 arrest for protein {protein.accession_number}")

                if pa.get(NORMALISED) and pa[NORMALISED].get(LOG2_MEAN):
                    self._import_protein_readings(
                        replicates_by_name,
                        sample_stages_by_name,
                        protein,
                        PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
                        pa[NORMALISED][LOG2_MEAN],
                    )
                else:
                    print(f"No normalised log2 mean for protein {protein.accession_number}")
                
                if pa.get(NORMALISED) and pa[NORMALISED].get(MIN_MAX):
                    self._import_protein_readings(
                        replicates_by_name,
                        sample_stages_by_name,
                        protein,
                        PROTEIN_ABUNDANCES_NORMALISED_MIN_MAX,
                        pa[NORMALISED][MIN_MAX],
                    )
                else:
                    print(f"No normalised min max for protein {protein.accession_number}")

                if pa.get(NORMALISED) and pa[NORMALISED].get(ZERO_MAX):
                    self._import_protein_readings(
                        replicates_by_name,
                        sample_stages_by_name,
                        protein,
                        PROTEIN_ABUNDANCES_NORMALISED_ZERO_MAX,
                        pa[NORMALISED][ZERO_MAX],
                    )
                else:
                    print(f"No normalised zero max for protein {protein.accession_number}")

                if pa.get(IMPUTED):
                    self._import_protein_readings(
                        replicates_by_name,
                        sample_stages_by_name,
                        protein,
                        PROTEIN_ABUNDANCES_IMPUTED,
                        pa[IMPUTED],
                    )
                else:
                    print(f"No imputed values for protein {protein.accession_number}")

                if pm.get(LOG2_MEAN):
                    self._import_metrics(
                        protein,
                        PROTEIN_ABUNDANCES_NORMALISED_LOG2_MEAN,
                        pm[LOG2_MEAN],
                    )
                else:
                    print(f"No log2 mean metrics {protein.accession_number}")



    def _import_metrics(
        self,
        protein,
        statistic_type_name,
        obj,
    ):
        _, stat_prot_raw = self._fetch_stats_type_and_stats(statistic_type_name, protein=protein)

        stat_prot_raw.metrics = self._convert_decimals(obj)
        stat_prot_raw.save()




    def _import_protein_readings(
        self,
        replicates_by_name,
        sample_stages_by_name,
        protein,
        statistic_type_name,
        obj
    ):
        _, stat_prot_raw = self._fetch_stats_type_and_stats(statistic_type_name, protein=protein)

        for replicate_name, readings in obj.items():
            for sample_stage_name, reading in readings.items():
                # Strip out 'status' string for imputed values
                if statistic_type_name == PROTEIN_ABUNDANCES_IMPUTED:
                    reading = reading['value']

                replicate = replicates_by_name[replicate_name]
                sample_stage = sample_stages_by_name[sample_stage_name]

                Abundance.objects.create(
                    statistic=stat_prot_raw,
                    replicate=replicate,
                    sample_stage=sample_stage,
                    reading=reading
                )



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

    def _convert_decimals(self, obj):
        if isinstance(obj, dict):
            return {k: self._convert_decimals(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_decimals(i) for i in obj]
        elif isinstance(obj, Decimal):
            return float(obj)
        return obj