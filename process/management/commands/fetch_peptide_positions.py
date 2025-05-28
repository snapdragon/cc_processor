import logging
from time import sleep

import pandas as pd
import requests
from django.core.management.base import BaseCommand

from process.models import PeptideStartPosition, Project

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "fetch peptide sequence start positions from uniprot"

    def handle(self, *args, **options):
        project_name = "SL"

        logger.info("Fetching peptide positions")

        project = Project.objects.get(name=project_name)

        file_path = f"data/{project.phosphoproteome_file}"

        excel_file = pd.ExcelFile(file_path)

        # Sheet two is the phosphoproteome sheet, apparently
        df = pd.read_excel(excel_file, sheet_name=1)

        row_no = 0

        for index, row in df.iterrows():
            accession_number = row[project.proteome_file_accession_number_column_name]

            if accession_number.startswith("CON_"):
                print(f"Skipping contaminant {accession_number}")
                continue

            if not row_no % 1000:
                logger.info(f"Adding phospho row {row_no} {accession_number}")

            row_no += 1

            peptide = row["Stripped.Sequence"]

            try:
                PeptideStartPosition.objects.get(
                    accession_number=accession_number, peptide=peptide
                )
            except Exception:
                print(f"Fetching position for {accession_number} {peptide}")

                start_position = get_peptide_start(accession_number, peptide)

                if start_position is None:
                    print(f"No start position found for {accession_number} {peptide}")
                    continue

                PeptideStartPosition.objects.create(
                    accession_number=accession_number,
                    peptide=peptide,
                    start_position=start_position,
                )

                sleep(0.01)


def get_peptide_start(uniprot_acc, peptide_seq):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError("Invalid UniProt accession or API error.")

    fasta = response.text
    protein_seq = "".join(fasta.split("\n")[1:])

    start_pos = protein_seq.find(peptide_seq)
    if start_pos == -1:
        return None
    return start_pos + 1
