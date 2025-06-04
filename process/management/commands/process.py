import logging
import math
import statistics

from django.core.management.base import BaseCommand, CommandError
from django.db.models.query import QuerySet

from process.models import ColumnName, Project, ProteinReading, Replicate

# TODO - make this configurable by flag
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "Processes all proteins for a given project"

    def add_arguments(self, parser):
        parser.add_argument(
            "--project",
            required=True,
            help="The name of the project to process",
        )
        parser.add_argument(
            "--with-bugs",
            help="Run with the original ICR bugs",
            action="store_true",
        )

    def handle(self, *args, **options):
        project_name = options["project"]
        with_bugs = options["with_bugs"]

        if with_bugs and project_name != "ICR":
            raise CommandError("Only an ICR project can run --with-bugs")

        logger.info(f"Processing for project {project_name}, with bugs {with_bugs}")

        project = Project.objects.get(name=project_name)
        replicates = Replicate.objects.filter(project=project)
        column_names = ColumnName.objects.filter(replicate__project=project)
        protein_readings = ProteinReading.objects.filter(
            column_name__replicate__project=project
        )

        self._process(project, replicates, protein_readings, column_names, with_bugs)

    def _process(
        self, project, replicates, protein_readings, column_names, with_bugs: bool
    ):
        """
        Does all the required calculations. The steps are:

        TODO - put the descriptions below in the function comments, not here, along with
            data structures.

        1) Calculate the median for each replicate for each stage.
            In other words, get all the values for each column and find the median.

        2) Calculate the mean for each protein for each stage across replicates.
            So for proteins 'ABCD', 'EFGH' with replicates 'One', 'Two' and stages '1h', '2h',
            it will get the abundances for 'ABCD 1h' for each of the replicates, then take
            the mean.

        3) Normalise all abundances. This is done by dividing each abundance by the median
            for its column, then averaging them across replicates.
        """

        # TODO - make each call flaggable
        # TODO - rename 'medians' to something more informative
        # TODO - does this need to be by replicate? Why not just all columns at once?
        medians = self._all_replicates(
            func=self._calc_replicate_column_medians,
            replicates=replicates,
            protein_readings=protein_readings,
            column_names=column_names,
        )

        if with_bugs:
            # Make both the same, it's simpler than deleting number one and changing
            #   code elsewhere
            replicate_1 = None
            replicate_2 = None

            for replicate in medians.keys():
                # TODO - make these constants
                if replicate.name == "One":
                    replicate_1 = replicate
                elif replicate.name == "Two":
                    replicate_2 = replicate

            r2_medians = {}

            for cn in medians[replicate_2].keys():
                r2_medians[cn.sample_stage.name] = medians[replicate_2][cn]

            for cn in medians[replicate_1].keys():
                medians[replicate_1][cn] = r2_medians[cn.sample_stage.name]

        print(f"Medians for project {project.name}")
        for replicate in medians.keys():
            print(f"Replicate: {replicate.name}")

            for stage in medians[replicate].keys():
                print(f"    {stage.name}: {medians[replicate][stage]}")

        # TODO - is this used for anything?
        # means_across_replicates_by_stage = (
        #     self._calculate_means_across_replicates_by_stage(
        #         protein_readings, with_bugs
        #     )
        # )

        # # These are just here for now to stop the pre commit hooks complaining
        # assert means_across_replicates_by_stage == means_across_replicates_by_stage

        normalised_protein_readings = self._calculate_first_level_normalisation(
            protein_readings, medians
        )

        logger.info(
            f"Number of normalised readings: f{len(normalised_protein_readings)}"
        )

        # TODO - check whether all means calculations need with-bugs
        normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                normalised_protein_readings, with_bugs
            )
        )

        assert (
            normalised_means_across_replicates_by_stage
            == normalised_means_across_replicates_by_stage
        )

        arrest_log2_normalised_protein_readings = (
            self._calculate_arrest_log2_normalisation(
                normalised_protein_readings, project
            )
        )

        logger.info(
            f"Number of arrest log2 normalised readings: f{len(arrest_log2_normalised_protein_readings)}"
        )

        arrest_log2_normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                arrest_log2_normalised_protein_readings, with_bugs
            )
        )

        assert (
            arrest_log2_normalised_means_across_replicates_by_stage
            == arrest_log2_normalised_means_across_replicates_by_stage
        )

        relative_log2_normalised_protein_readings = (
            self._calculate_relative_log2_normalisation(
                normalised_protein_readings, project
            )
        )

        logger.info(
            f"Number of relative log2 normalised readings: f{len(relative_log2_normalised_protein_readings)}"
        )

        # relative_log2_normalised_means_across_replicates_by_stage = (
        #     self._calculate_means_across_replicates_by_stage(
        #         relative_log2_normalised_protein_readings, with_bugs
        #     )
        # )

    def _all_replicates(self, *args, **kwargs):
        """
        Calls the passed function for each replicate for the project.
        """
        results = {}

        # Remove the passed function and replicates as they're not needed by the passed function
        call_kwargs = kwargs.copy()
        func = call_kwargs.pop("func")
        replicates = call_kwargs.pop("replicates")

        for replicate in replicates:
            call_kwargs["replicate"] = replicate
            results[replicate] = func(**call_kwargs)

        return results

    def _calculate_relative_log2_normalisation(
        self, normalised_protein_readings: QuerySet[ProteinReading], project: Project
    ):
        logger.info("relative log2 normalising abundances")

        log2_normalised_protein_readings: list = []

        for protein_reading in normalised_protein_readings:
            if (
                protein_reading.protein.accession_number == "Q09666"
                and protein_reading.reading is not None
            ):
                log2_reading = math.log2(protein_reading.reading)

                print(
                    f"log2 norm for Q09666 {protein_reading.column_name.sample_stage.name}: {log2_reading}"
                )

            # log2_normalised_protein_readings.append()

        return log2_normalised_protein_readings

    def _calculate_arrest_log2_normalisation(
        self, normalised_protein_readings: QuerySet[ProteinReading], project: Project
    ):
        logger.info("log2 arrest normalising abundances")

        log2_normalised_protein_readings = []

        # TODO - what should the stage name be?
        ARRESTING_AGENT = "Nocodozole"

        # TODO - this is a hack, maybe add the field to the Project model?
        if project.name == "ICR":
            ARRESTING_AGENT = "Palbo"

        protein_arrest_values: dict[Replicate, dict] = {}

        # Get all arrest values by replicate by protein
        # TODO - this is super inefficient
        for protein_reading in normalised_protein_readings:
            if protein_reading.column_name.sample_stage.name == ARRESTING_AGENT:
                if not protein_arrest_values.get(protein_reading.column_name.replicate):
                    protein_arrest_values[protein_reading.column_name.replicate] = {}

                protein_arrest_values[protein_reading.column_name.replicate][
                    protein_reading.protein
                ] = protein_reading.reading

        for protein_reading in normalised_protein_readings:
            log2_reading = None

            if (
                protein_reading.reading is not None
                and protein_arrest_values[protein_reading.column_name.replicate][
                    protein_reading.protein
                ]
                is not None
            ):
                log2_reading = self._round(
                    math.log2(
                        protein_reading.reading
                        / protein_arrest_values[protein_reading.column_name.replicate][
                            protein_reading.protein
                        ]
                    )
                )

            log2_normalised_protein_readings.append(
                ProteinReading(
                    protein=protein_reading.protein,
                    column_name=protein_reading.column_name,
                    reading=log2_reading,
                )
            )

        return log2_normalised_protein_readings

    def _calculate_first_level_normalisation(
        self, protein_readings: QuerySet[ProteinReading], medians
    ):
        logger.info("Normalising abundances")

        normalised_protein_readings = []

        protein_no = 0

        for protein_reading in protein_readings:
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Normalising for {protein_no}, {protein_reading.protein.accession_number} {protein_reading.column_name.sample_stage.name}",
            )

            reading = protein_reading.reading

            normalised_reading = None

            if reading:
                median = medians[protein_reading.column_name.replicate][
                    protein_reading.column_name
                ]

                # TODO - need to round this?
                # TODO - might there be undefined medians? If so what to do?
                #   Maybe throw an exception in the median calculation func if one is encountered?
                #   They probably shouldn't happen.
                normalised_reading = reading / median

            # TODO - not a QuerySet
            # TODO - inefficient
            normalised_protein_readings.append(
                ProteinReading(
                    protein=protein_reading.protein,
                    column_name=protein_reading.column_name,
                    reading=normalised_reading,
                )
            )

        return normalised_protein_readings

    def _calculate_means_across_replicates_by_stage(
        self, protein_readings: QuerySet[ProteinReading], with_bugs: bool
    ):
        logger.info("Calculating mean across replicates by stage for each protein")

        means: dict = {}

        protein_no = 0

        for protein_reading in protein_readings:
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Calculating mean for {protein_no}, {protein_reading.protein.accession_number} {protein_reading.column_name.sample_stage.name}",
            )

            protein = protein_reading.protein

            if not means.get(protein):
                means[protein] = {}

            stage = protein_reading.column_name.sample_stage

            if not means[protein].get(stage):
                means[protein][stage] = []

            if with_bugs:
                # We throw away the second reading
                if len(means[protein][stage]) == 1:
                    continue

            if protein_reading.reading is not None:
                means[protein][stage].append(protein_reading.reading)

        protein_stage_no = 0

        for protein in means.keys():
            for stage in means[protein]:
                protein_stage_no += 1

                abundances = means[protein][stage]

                mean = None

                if len(abundances):
                    mean = sum(abundances) / len(abundances)
                    means[protein][stage] = self._round(mean)
                else:
                    # TODO - is this the right thing to do?
                    means[protein][stage] = None

                protein_stage_no += 1
                self._count_logger(
                    protein_stage_no,
                    10000,
                    f"Mean for {protein_stage_no}, {protein} {stage}: {mean}",
                )

        return means

    def _calc_replicate_column_medians(
        self,
        replicate: Replicate,
        protein_readings: QuerySet[ProteinReading],
        column_names: QuerySet[ColumnName],
    ):
        logger.info("Calculating column medians by replicate")

        column_medians = {}

        column_names_by_replicate = column_names.filter(replicate__name=replicate.name)

        for column_name in column_names_by_replicate:
            readings = []

            protein_readings_by_column = protein_readings.filter(
                column_name=column_name
            )

            for protein_reading in protein_readings_by_column:
                # TODO - what to do about None values?
                if protein_reading.reading is not None:
                    readings.append(protein_reading.reading)

            median = statistics.median(readings)

            column_medians[column_name] = median

        return column_medians

    def _count_logger(self, i: int, step: int, output: str):
        if i % step == 0:
            logger.info(output)

    def _round(self, value):
        return round(value, 4)
