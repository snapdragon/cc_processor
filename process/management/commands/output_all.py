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
    help = "Output a protein's combine results to a file."

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to output for.",
        )
        parser.add_argument(
            "--with-bugs",
            help="Run with the original ICR bugs",
            action="store_true",
        )

    def handle(self, *args, **options):
        project_name = options["project"]
        with_bugs = options["with_bugs"]

        if not project_name:
            raise Exception("Please provide a project name")

        project = Project.objects.get(name=project_name)

        logger.info(f"Outputting all proteins for {project_name} with bugs {with_bugs}")

        run_results = RunResult.objects.get(run__project=project, run__with_bugs=with_bugs)

        all_results = {}

        for rr in run_results:
            all_results[rr.protein.accession_number] = rr.protein_phospho_results

        file_path = f"output/{project_name}_with_bugs_{with_bugs}_all.json"

        with open(file_path, "w") as outfile:
            json.dump(all_results, outfile, indent=4)
