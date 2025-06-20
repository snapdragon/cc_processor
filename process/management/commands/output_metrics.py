import logging
import json

import pandas as pd
from django.core.management.base import BaseCommand

from process.models import Project, Protein, RunResult

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "Output a protein's metrics."

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to output for.",
        )
        parser.add_argument(
            "--accession-number",
            required=True,
            help="The accession number of the protein to output.",
        )
        parser.add_argument(
            "--with-bugs",
            help="Run with the original ICR bugs",
            action="store_true",
        )

    def handle(self, *args, **options):
        project_name = options["project"]
        accession_number = options["accession_number"]
        with_bugs = options["with_bugs"]

        if not project_name:
            raise Exception("Please provide a project name")

        if not accession_number:
            raise Exception("Please provide an accession number")

        project = Project.objects.get(name=project_name)
        protein = Protein.objects.get(accession_number=accession_number, project=project)

        logger.info(f"Outputting protein {protein} for {project_name}")

        run_result = RunResult.objects.get(protein=protein, run__project=project, run__with_bugs=with_bugs)

        file_path = f"output/{project_name}_{protein}_with_bugs_{with_bugs}_metrics.json"

        with open(file_path, "w") as outfile:
            json.dump(run_result.protein_phospho_result['metrics'], outfile, indent=4)
