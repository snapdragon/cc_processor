import logging

import pandas as pd
from django.core.management.base import BaseCommand

from process.models import ColumnName, Project, Protein, ProteinReading

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

# TODO - this is a duplicate, move somewhere common
PROJECT_NAME = "Soliman Labs"


class Command(BaseCommand):
    help = "import a spreadsheet"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            default=PROJECT_NAME,
            help="Project name",
        )

    def handle(self, *args, **options):
        project_name = options["project"]

        logger.info("Importing spreadsheet for {project_name}")

        project = Project.objects.get(name=project_name)

        file_path = f"data/{project.proteome_file}"

        excel_file = pd.ExcelFile(file_path)

        df = pd.read_excel(excel_file, sheet_name=0)

        column_names = ColumnName.objects.all()

        cns_by_name = {}

        # TODO - take this out, or modify it to the project, when adding more projects
        Protein.objects.all().delete()
        ProteinReading.objects.all().delete()

        for cn in column_names:
            cns_by_name[cn.name] = cn

        row_no = 0

        for index, row in df.iterrows():
            row_no += 1

            print(f"Adding row {row_no}")

            accession_number = row[project.proteome_file_accession_number_column_name]

            protein = Protein.objects.create(accession_number=accession_number)

            for col in df.columns:
                if cn := cns_by_name.get(col):
                    reading = row[col]

                    if reading != reading:
                        reading = None

                    ProteinReading.objects.create(
                        column_name=cn, reading=reading, protein=protein
                    )
