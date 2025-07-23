import logging
from time import sleep
import json

from django.core.management.base import BaseCommand

from process.models import Project, Protein, Phospho

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "fetch peptide sequence start positions from uniprot"

    def handle(self, *args, **options):
        logger.info("Export supplementary data for Soliman Labs")

        project = Project.objects.get(name="SL")

        proteins = Protein.objects.filter(project=project, relevant=True)

        output_file_path = 'supplementary_data/relevant_proteins.txt'

        data = []

        for protein in proteins:
            data.append(protein.accession_number)

        # Write to file
        with open(output_file_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)

        phosphos = Phospho.objects.filter(protein__project=project, relevant=True)

        output_file_path = 'supplementary_data/relevant_phospho_proteins.txt'

        data = []

        for phospho in phosphos:
            data.append(phospho.protein.accession_number)

        # Write to file
        with open(output_file_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)
