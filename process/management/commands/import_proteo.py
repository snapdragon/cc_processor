import logging

import pandas as pd
from django.core.management.base import BaseCommand

from process.models import ColumnName, Project, Protein, ProteinReading

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "import a spreadsheet"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to import",
        )

    def handle(self, *args, **options):
        project_name = options["project"]

        if not project_name:
            raise Exception("Please provide a project name")

        logger.info(f"Importing spreadsheet for {project_name}")

        project = Project.objects.get(name=project_name)

        file_path = f"data/{project.proteome_file}"

        excel_file = pd.ExcelFile(file_path)

        df = pd.read_excel(excel_file, sheet_name=0)

        # TODO - untested, may fail
        column_names = ColumnName.objects.filter(replicate__project=project)
        cns_by_name = {}
        for cn in column_names:
            cns_by_name[cn.name] = cn

        # TODO - take this out, or modify it to the project, when adding more projects
        Protein.objects.filter(project=project).delete()
        ProteinReading.objects.filter(protein__project=project).delete()


        row_no = 0

        for index, row in df.iterrows():
            row_no += 1

            is_contaminant = False

            # TODO - this is a hack for ICR. Maybe make it a DB config?
            if contaminant := row.get("Contaminant"):
                if contaminant == "Yes" or contaminant == "TRUE":
                    is_contaminant = True

            print(f"Adding row {row_no}")

            accession_number = row[project.proteome_file_accession_number_column_name]

            protein = Protein.objects.create(
                project=project, accession_number=accession_number, is_contaminant=is_contaminant
            )

            # We don't want the readings for contaminants
            if is_contaminant:
                continue

            for col in df.columns:
                if cn := cns_by_name.get(col):
                    reading = row[col]

                    if reading != reading:
                        reading = None

                    ProteinReading.objects.create(
                        column_name=cn, reading=reading, protein=protein
                    )

        print(f"Total rows: {row_no}")
