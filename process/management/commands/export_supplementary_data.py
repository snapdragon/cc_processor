import logging
from time import sleep
import json

from django.core.management.base import BaseCommand

from process.models import Project, Protein, Phospho, Abundance

from process.constants import ABUNDANCES_RAW

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


        # Relevant proteins and phosphos

        proteins = Protein.objects.filter(project=project, relevant=True)

        output_file_path = 'supplementary_data/relevant_proteins.txt'

        data = []

        for protein in proteins:
            data.append(protein.accession_number)

        with open(output_file_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)

        phosphos = Phospho.objects.filter(protein__project=project, relevant=True)

        output_file_path = 'supplementary_data/relevant_phospho_proteins.txt'

        data = []

        for phospho in phosphos:
            data.append(phospho.protein.accession_number)

        with open(output_file_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)



        # Output GO enrichment lists

        output_file_path = 'supplementary_data/protein_GO_enrichment_list.txt'

        with open(output_file_path, 'w', encoding='utf-8') as f:
            json.dump(project.protein_go_list, f, ensure_ascii=False, indent=4)

        output_file_path = 'supplementary_data/phospho_protein_GO_enrichment_list.txt'

        with open(output_file_path, 'w', encoding='utf-8') as f:
            json.dump(project.phospho_protein_go_list, f, ensure_ascii=False, indent=4)



        # Output percentages of blank cells for proteins and phosphos

        protein_rows = Protein.objects.filter(project=project).count()
        phospho_rows = Phospho.objects.filter(protein__project=project).count()

        raw_protein_abundances_num = Abundance.objects.filter(
            statistic__protein__project = project,
            statistic__statistic_type__name = ABUNDANCES_RAW,
            replicate__mean = False
        ).count()

        raw_phospho_abundances_num = Abundance.objects.filter(
            statistic__phospho__protein__project = project,
            statistic__statistic_type__name = ABUNDANCES_RAW,
            replicate__mean = False
        ).count()

        protein_blank_percentage = round((raw_protein_abundances_num / (protein_rows * 30)) * 100, 2)
        phospho_blank_percentage = round((raw_phospho_abundances_num / (phospho_rows * 30)) * 100, 2)

        print(f"Total populated proteins: {raw_protein_abundances_num} of {protein_rows * 30}: {protein_blank_percentage}% populated")
        print(f"Total populated phosphos: {raw_phospho_abundances_num} of {phospho_rows * 30}: {phospho_blank_percentage}% populated")
