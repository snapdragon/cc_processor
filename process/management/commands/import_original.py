import logging

import pandas as pd
from django.core.management.base import BaseCommand
import ijson

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
    PROTEIN_ABUNDANCES_RAW
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
                    print(f"Processing {gene_name}")

                protein = Protein.objects.create(
                    project=project, accession_number=gene_name, is_contaminant=False
                )

                _, stat_prot_raw = self._fetch_stats_type_and_stats(PROTEIN_ABUNDANCES_RAW, protein=protein)

                pa = gene_data[PROTEIN_ABUNDANCES]

                for replicate_name, readings in pa[RAW].items():
                    for sample_stage_name, reading in readings.items():
                        replicate = replicates_by_name[replicate_name]
                        sample_stage = sample_stages_by_name[sample_stage_name]

                        Abundance.objects.create(
                            statistic=stat_prot_raw,
                            replicate=replicate,
                            sample_stage=sample_stage,
                            reading=reading
                        )

                # if protein.accession_number == 'Q93075':
                #     print(Abundance.objects.filter(statistic=stat_prot_raw))




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
