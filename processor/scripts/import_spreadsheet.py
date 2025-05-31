import os
import sys

import django
import pandas as pd

# Add project root to PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# Set Django settings module
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "processor.settings")

django.setup()

from process.models import (ColumnName, Project, Protein,  # noqa: E402
                            ProteinReading)

project = Project.objects.get(name="Soliman Labs")

file_path = f"data/{project.proteome_file}"

excel_file = pd.ExcelFile(file_path)

df = pd.read_excel(excel_file, sheet_name=0)

column_names = ColumnName.objects.all()

cns_by_name = {}

for cn in column_names:
    cns_by_name[cn.name] = cn

# TODO - take this out, or modify it to the project, when adding more projects
Protein.objects.all().delete()
ProteinReading.objects.all().delete()

for index, row in df.iterrows():
    accession_number = row["Protein.Group"]

    protein = Protein.objects.create(accession_number=accession_number)

    for col in df.columns:
        if cn := cns_by_name.get(col):
            ProteinReading.objects.create(column_name=cn, reading=row[col])
