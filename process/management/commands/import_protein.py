import logging
import re

import pandas as pd
from django.core.management.base import BaseCommand

from process.models import ColumnName, Project, Protein, Abundance, Phospho, StatisticType, Statistic

from process.constants import (
    ABUNDANCES_RAW
)

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

        # Change script name

        # TODO - output available project names
        if not project_name:
            raise Exception(f"Invalud project name {project_name}")

        logger.info(f"Importing spreadsheet for {project_name}")

        project = Project.objects.get(name=project_name)
        stats_type_rp = StatisticType.objects.get(name=ABUNDANCES_RAW)

        file_path = f"data/{project.proteome_file}"

        excel_file = pd.ExcelFile(file_path)

        df = pd.read_excel(excel_file, sheet_name=0)

        column_names = ColumnName.objects.filter(replicate__project=project)
        cns_by_name = {}
        for cn in column_names:
            cns_by_name[cn.name] = cn

        Abundance.objects.filter(statistic__protein__project=project).delete()
        Statistic.objects.filter(protein__project=project).delete()
        Protein.objects.filter(project=project).delete()
        Phospho.objects.filter(protein__project=project).delete()

        row_no = 0

        for index, row in df.iterrows():
            if row_no % 1000 == 0:
                print(f"Adding row {row_no}")

            row_no += 1

            is_contaminant = False

            # TODO - this is a hack for ICR. Maybe make it a DB config?
            if contaminant := row.get("Contaminant"):
                if contaminant == "Yes" or contaminant == "TRUE":
                    is_contaminant = True

            accession_number = row[project.proteome_file_accession_number_column_name]

            # SL contaminants start with 'CON_'. I think.
            #   They can't be looked up in uniprot anyway.
            if accession_number.startswith("CON_"):
                is_contaminant = True

            # We don't want the readings for contaminants
            if is_contaminant:
                print(f"Skipping contaminant {accession_number}")
                continue

            if protein := Protein.objects.filter(project=project, accession_number=accession_number).first():
                print(f"DUPLICATE ACCESSION NUMBER {accession_number}")

                # To ensure the two readings don't get mashed up if they
                #   have any missing values
                Abundance.objects.filter(statistic__protein=protein).delete()

                statistic = Statistic.objects.get(
                    statistic_type=stats_type_rp, protein=protein
                )
            else:
                protein = Protein.objects.create(
                    project=project, accession_number=accession_number, is_contaminant=is_contaminant
                )

                statistic = Statistic.objects.create(
                    statistic_type=stats_type_rp, protein=protein
                )

            if project_name == "SL":
                for col in df.columns:
                    col_short = re.search(r'IITI_\d{3}', col)

                    if col_short is not None:
                        if cn := cns_by_name.get(col_short.group()):
                            reading = row[col]

                            if reading != reading:
                                continue

                            Abundance.objects.create(
                                statistic=statistic, replicate=cn.replicate, sample_stage=cn.sample_stage, reading=reading
                            )
            else:
                for col in df.columns:
                    if cn := cns_by_name.get(col):
                        reading = row[col]

                        if reading != reading:
                            continue

                        Abundance.objects.create(
                            statistic=statistic, replicate=cn.replicate, sample_stage=cn.sample_stage, reading=reading
                        )


        print(f"Total rows: {row_no}")
