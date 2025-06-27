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

def compare_json(obj1, obj2, path=""):
    differences = []

    if isinstance(obj1, dict) and isinstance(obj2, dict):
        all_keys = set(obj1.keys()).union(obj2.keys())
        for key in all_keys:
            new_path = f"{path}.{key}" if path else key
            val1 = obj1.get(key, "__MISSING__")
            val2 = obj2.get(key, "__MISSING__")
            differences += compare_json(val1, val2, new_path)

    elif isinstance(obj1, list) and isinstance(obj2, list):
        for i, (item1, item2) in enumerate(zip(obj1, obj2)):
            new_path = f"{path}[{i}]"
            differences += compare_json(item1, item2, new_path)
        if len(obj1) != len(obj2):
            differences.append(f"{path}: List lengths differ ({len(obj1)} != {len(obj2)})")

    else:
        if obj1 != obj2:
            differences.append(f"{path}: {obj1} != {obj2}")

    return differences

class Command(BaseCommand):
    help = "Compare output files from ICR and Process"

    def add_arguments(self, parser):
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
        accession_number = options["accession_number"]
        with_bugs = options["with_bugs"]

        if not accession_number:
            raise Exception("Please provide an accession number")

        logger.info(f"Comparing files for  {accession_number} with_bugs {with_bugs}")

        process_file_path = f"output/ICR_{accession_number}_with_bugs_{with_bugs}_combined_result.json"
        ICR_file_path = f"output_ICR/TimeCourse_{accession_number}_info.json"

        # with open(process_file_path) as file:
        #     process_result = json.load(process_file_path)

        # with open(ICR_file_path) as file:
        #     ICR_result = json.load(process_file_path)

        with open(process_file_path) as f1, open(ICR_file_path) as f2:
            process_result = json.load(f1)
            ICR_result = json.load(f2)

        diffs = compare_json(process_result, ICR_result)

        for diff in diffs:
            print(diff)
