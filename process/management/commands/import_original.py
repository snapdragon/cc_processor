import logging

import pandas as pd
from django.core.management.base import BaseCommand
import ijson

from process.models import Project, Protein, Abundance, Phospho, StatisticType, Statistic

from process.constants import (
    PROTEIN_READINGS
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
        stats_type_rp = StatisticType.objects.get(name=PROTEIN_READINGS)

        file_path = f"data/ICR/TimeCourse_Full_info_full_indented.json"

        Abundance.objects.filter(statistic__protein__project=project).delete()
        Statistic.objects.filter(protein__project=project).delete()
        Protein.objects.filter(project=project).delete()
        Phospho.objects.filter(protein__project=project).delete()

        with open(file_path, 'r') as f:
            for gene_name, gene_data in ijson.kvitems(f, ''):
                print(f"Processing {gene_name}")
                
        # for index, row in df.iterrows():
        #     if row_no % 1000 == 0:
        #         print(f"Adding row {row_no}")

        #     row_no += 1

        #     is_contaminant = False

        #     # TODO - this is a hack for ICR. Maybe make it a DB config?
        #     if contaminant := row.get("Contaminant"):
        #         if contaminant == "Yes" or contaminant == "TRUE":
        #             is_contaminant = True

        #     accession_number = row[project.proteome_file_accession_number_column_name]

        #     if protein := Protein.objects.filter(project=project, accession_number=accession_number).first():
        #         print(f"DUPLICATE ACCESSION NUMBER {accession_number}")

        #         # To ensure the two readings don't get mashed up if they
        #         #   have any missing values
        #         Abundance.objects.filter(statistic__protein=protein).delete()

        #         statistic = Statistic.objects.get(
        #             statistic_type=stats_type_rp, protein=protein
        #         )
        #     else:
        #         protein = Protein.objects.create(
        #             project=project, accession_number=accession_number, is_contaminant=is_contaminant
        #         )

        #         statistic = Statistic.objects.create(
        #             statistic_type=stats_type_rp, protein=protein
        #         )

        #     # We don't want the readings for contaminants
        #     if is_contaminant:
        #         continue

        #     for col in df.columns:
        #         if cn := cns_by_name.get(col):
        #             reading = row[col]

        #             if reading != reading:
        #                 continue

        #             Abundance.objects.create(
        #                 statistic=statistic, replicate=cn.replicate, sample_stage=cn.sample_stage, reading=reading
        #             )

        # print(f"Total rows: {row_no}")
