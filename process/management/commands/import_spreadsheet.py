import pandas as pd
from django.core.management.base import BaseCommand

from process.models import ColumnName, Project, Protein, ProteinReading


class Command(BaseCommand):
    help = "import a spreadsheet"

    def handle(self, *args, **kwargs):
        project = Project.objects.get(name="Soliman Labs")

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

            accession_number = row["Protein.Group"]

            protein = Protein.objects.create(accession_number=accession_number)

            for col in df.columns:
                if cn := cns_by_name.get(col):
                    reading = row[col]

                    if reading != reading:
                        reading = None

                    ProteinReading.objects.create(
                        column_name=cn, reading=reading, protein=protein
                    )
