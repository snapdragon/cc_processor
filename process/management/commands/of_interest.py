import json

from django.core.management.base import BaseCommand

from process.constants import (
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ANOVA,
    CURVE_FOLD_CHANGE,
    Q_VALUE,
)
from process.models import Project, Statistic

OUTPUT_DIR = "output/lists"


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=False,
            help="The name of the project to import",
        )
        parser.add_argument(
            "--accession-number",
            required=False,
            help="The accession number of the protein to compare.",
        )

    def handle(self, *args, **options):
        print("Listing proteins and phosphos of interest")

        project_name = options["project"]

        if project_name:
            project = Project.objects.get(name=project_name)

            self._run_project(project)
            return

        projects = Project.objects.all()

        all_matches = {}

        for project in projects:
            all_matches[project.name] = self._run_project(project)

        for pr1 in all_matches:
            for pr2 in all_matches:
                if pr1 == pr2:
                    continue

                print(f"Matching {pr1} to {pr2}")

                am1 = all_matches[pr1][0]
                am2 = all_matches[pr2][0]

                res = [value for value in am1 if value in am2]
                print(f"Protein overlap {len(res)}")

                am1 = all_matches[pr1][1]
                am2 = all_matches[pr2][1]

                res = [value for value in am1 if value in am2]
                print(f"Phospho overlap {len(res)}")

                print("")

    def _run_project(self, project):
        if project.name == "ICR":
            print(f"Project: {project.name} with_bugs {project.with_bugs}")
        else:
            print(f"Project: {project.name}")

        # if accession_number:
        #     statistics = Statistic.objects.filter(
        #         statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
        #         protein__project=project,
        #         protein__accession_number=accession_number
        #     )
        # else:
        statistics = Statistic.objects.filter(
            statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
            protein__project=project,
        )

        protein_matches = self._run_stats(statistics, project, 0.05)

        # if accession_number:
        #     statistics = Statistic.objects.filter(
        #         statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
        #         phospho__protein__project=project,
        #         phospho__protein__accession_number=accession_number
        #     )
        # else:
        statistics = Statistic.objects.filter(
            statistic_type__name=ABUNDANCES_NORMALISED_LOG2_MEAN,
            phospho__protein__project=project,
        )

        phospho_matches = self._run_stats(statistics, project, 0.02)

        print("")

        return protein_matches, phospho_matches

    def _run_stats(self, statistics, project, max_q):
        total = 0
        matches = set()
        num_matches = 0

        for statistic in statistics:
            total += 1

            if (
                statistic.metrics is None
                or statistic.metrics.get(ANOVA) is None
                or statistic.metrics[ANOVA].get(Q_VALUE) is None
                or statistic.metrics[ANOVA][Q_VALUE] >= max_q
            ):
                # print(f"q value {statistic.metrics[ANOVA][Q_VALUE]}")
                continue

            if (
                statistic.metrics.get(CURVE_FOLD_CHANGE) is None
                or statistic.metrics[CURVE_FOLD_CHANGE] < 1.2
            ):
                continue

            num_matches += 1

            if statistic.protein:
                matches.add(statistic.protein.accession_number)
            else:
                matches.add(statistic.phospho.protein.accession_number)

            # if statistic.protein:
            #     print(f"Match for {statistic.protein.accession_number}")
            # else:
            #     print(f"Match for {statistic.phospho.protein.accession_number}")

        if statistic.protein:
            output = "Total protein"
        else:
            output = "Total phospho"

        print(
            f"{output} matches {len(matches)} of {total} ({num_matches} not distinct)"
        )

        file = project.name

        if project.name == "ICR":
            file = f"{file}_with_bugs_{project.with_bugs}"

        if statistic.protein:
            file = f"{file}_protein"
        else:
            file = f"{file}_phospho"

        with open(f"{OUTPUT_DIR}/{file}.json", "w") as f:
            json.dump(list(matches), f, indent=4)

        return matches
