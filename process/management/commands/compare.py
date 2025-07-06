import logging
import json

import pandas as pd
from django.core.management.base import BaseCommand

from process.models import (
    Project,
    Protein,
    Abundance,
    Phospho,
    StatisticType,
    Statistic,
    Replicate,
    SampleStage
)

from process.constants import (
    PROTEIN_ABUNDANCES_RAW,
    RAW,
    PROTEIN_ABUNDANCES,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--accession-number",
            required=True,
            help="The accession number of the protein to output.",
        )

    def handle(self, *args, **options):
        accession_number = options["accession_number"]

        if not accession_number:
            raise Exception("Please provide an accession number")

        logger.info(f"Comparing original and new")

        project_original = Project.objects.get(name="ICR")
        project_process = Project.objects.get(name="Original")

        if not project_original.with_bugs:
            raise Exception("Only a with_bugs ICR project can be compared to the original. You'll have to run 'process' again with --with-bugs.")

        # replicates_original = Replicate.objects.filter(
        #     project=project_original
        # )
        # sample_stages_original = SampleStage.objects.filter(
        #     project=project_original
        # ).order_by('rank')

        # replicates_process = Replicate.objects.filter(
        #     project=project_process
        # )
        # sample_stages_process = SampleStage.objects.filter(
        #     project=project_process
        # ).order_by('rank')

        protein_original = Protein.objects.get(
            project=project_original,
            accession_number = accession_number
        )
        protein_process = Protein.objects.get(
            project=project_process,
            accession_number = accession_number
        )

        _, stat_prot_raw_original = self._fetch_stats_type_and_stats(
            PROTEIN_ABUNDANCES_RAW,
            protein=protein_original
        )
        _, stat_prot_raw_process = self._fetch_stats_type_and_stats(
            PROTEIN_ABUNDANCES_RAW,
            protein=protein_process
        )

        stat_prot_raw_abundances_original = Abundance.objects.filter(
            statistic = stat_prot_raw_original
        )
        stat_prot_raw_abundances_process = Abundance.objects.filter(
            statistic = stat_prot_raw_process
        )

        print(len(stat_prot_raw_abundances_original))
        print(len(stat_prot_raw_abundances_process))

        for ab_original in stat_prot_raw_abundances_original:
            ab_process = stat_prot_raw_abundances_process.filter(
                replicate__name = ab_original.replicate.name,
                sample_stage__name = ab_original.sample_stage.name
            ).first()

            if ab_original.reading != ab_process.reading:
                print(f"No match for {ab_original.replicate.name} for {ab_original.sample_stage.name}")
            else:
                print(f"Raw matched {protein_original.accession_number} reading {ab_original.reading}")


    # TODO - copied from import_original
    def _fetch_stats_type_and_stats(self, statistic_type_name, project = None, protein = None, phospho = None):
        statistic_type = StatisticType.objects.get(name=statistic_type_name)

        if project:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, project=project)
        elif protein:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, protein=protein)
        elif phospho:
            stat, _ = Statistic.objects.get_or_create(statistic_type=statistic_type, phospho=phospho)
        else:
            raise Exception(f"_clear_and_fetch_stats needs a project, protein or phospho")

        return statistic_type, stat
