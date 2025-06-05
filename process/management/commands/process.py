import logging
import math
import statistics

import numpy as np
from django.core.management.base import BaseCommand, CommandError
from django.db.models.query import QuerySet
from scipy.stats import moment
from sklearn.metrics import r2_score

from process.models import ColumnName, Project, Protein, ProteinReading, Replicate

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
        # TODO - this is solely for development. Take it out afterwards.
        Q09666 = Protein.objects.get(
            accession_number="Q09666", project__name=project.name
        )

        # TODO - make each call flaggable
        # TODO - rename 'medians' to something more informative
        # TODO - does this need to be by replicate? Why not just all columns at once?
        # TODO - is _all_replicates really useful?
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

        # N.B. normalised_protein_readings is not the same structure as protein_readings.
        #   protein_readings is just a list of ProteinReading objects. normalised_protein_readings
        #   is a dict with Protein object keys. The values is a dict of replicate name keys
        #   with a dict of sample stage names and abundances as keys.
        normalised_protein_readings = self._calculate_first_level_normalisation(
            protein_readings, medians
        )

        logger.info(
            f"Number of normalised readings: f{len(normalised_protein_readings.keys())}"
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

        relative_log2_normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                relative_log2_normalised_protein_readings, with_bugs
            )
        )

        print("+++++ RELATIVE LOG2 NORM ABUNDANCES AVERAGES")
        print(relative_log2_normalised_means_across_replicates_by_stage[Q09666])

        level_two_normalised_protein_readings = self._calculate_level_two_normalisation(
            relative_log2_normalised_protein_readings
        )

        logger.info(
            f"Number of level two normalised readings: f{len(level_two_normalised_protein_readings)}"
        )

        level_two_normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                level_two_normalised_protein_readings, with_bugs
            )
        )

        assert (
            level_two_normalised_means_across_replicates_by_stage
            == level_two_normalised_means_across_replicates_by_stage
        )

        min_max_normalised_protein_readings = self._calculate_level_two_normalisation(
            normalised_protein_readings, True
        )

        logger.info(
            f"min max normalised readings: f{len(min_max_normalised_protein_readings)}"
        )

        min_max_normalised_means_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                min_max_normalised_protein_readings, with_bugs
            )
        )

        assert (
            min_max_normalised_means_across_replicates_by_stage
            == min_max_normalised_means_across_replicates_by_stage
        )

        # print("++++++ min max normalised means across replicates")
        # print(min_max_normalised_means_across_replicates_by_stage[Q09666])

        imputed_protein_readings = self._impute(
            level_two_normalised_protein_readings, replicates, column_names
        )

        logger.info(f"imputed readings: f{len(imputed_protein_readings)}")

        imputed_across_replicates_by_stage = (
            self._calculate_means_across_replicates_by_stage(
                imputed_protein_readings, with_bugs, imputed=True
            )
        )

        # print("++++++ imputed averages")
        # print(imputed_across_replicates_by_stage[Q09666])
        assert imputed_across_replicates_by_stage == imputed_across_replicates_by_stage

        log2_mean_metrics = self._calculate_metrics(
            relative_log2_normalised_protein_readings,
            relative_log2_normalised_means_across_replicates_by_stage,
        )

        print("+++++ LOG2 MEAN METRICS")
        print(log2_mean_metrics[Q09666])

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

    def _calculate_metrics(
        self,
        protein_readings: QuerySet[ProteinReading],
        protein_readings_averages: QuerySet[Replicate],
    ):
        logger.info("Calculate metrics")

        protein_metrics = {}

        for protein in protein_readings.keys():
            abundances = []
            abundance_averages = protein_readings_averages[protein]
            # TODO - how does ICR cope with Nones but not this? Does it use imputed values?
            # abundance_averages_list = [val for val in abundance_averages.values() if val is not None]
            abundance_averages_list = list(abundance_averages.values())

            for replicate_name in protein_readings[protein]:
                for stage_name in protein_readings[protein][replicate_name]:
                    abundance = protein_readings[protein][replicate_name][stage_name]

                    if abundance is not None:
                        # TODO - why does this add all reps for standard deviation calculation?
                        abundances.append(abundance)

            std_dev = statistics.stdev(abundances)

            try:
                _, residuals, *_ = np.polyfit(
                    range(len(abundance_averages)),
                    abundance_averages_list,
                    2,
                    full=True,
                )

                if len(residuals) == 0:
                    # TODO - why 5?
                    # TODO - 5 for ICR, but what about others?
                    # eg Q9HBL0 {'G2_2': 0.4496, 'G2/M_1': 0.7425, 'M/Early G1': 1.0}
                    residual = 5
                else:
                    residual = residuals[0]

                # TODO - find out what this all means
                r_squared = self._polyfit(
                    range(0, len(abundance_averages)), abundance_averages_list, 2
                )
                # TODO - find out what this all means
                max_fold_change = max(abundance_averages.values()) - min(
                    # TODO - change to abundance_averages_list?
                    abundance_averages.values()
                )

                # TODO - what does all this mean?
                base_metrics = {
                    "variance": moment(abundance_averages_list, moment=2),
                    "skewness": moment(abundance_averages_list, moment=3),
                    "kurtosis": moment(abundance_averages_list, moment=4),
                    "peak": max(abundance_averages, key=abundance_averages.get),
                    "max_fold_change": max_fold_change,
                    "residuals": residual,
                    "R_squared": r_squared,
                }

                # If we have info for the protein in at least 2 replicates
                if len(protein_readings[protein]) >= 2:
                    curve_fold_change, curve_peak = self._calcCurveFoldChange(
                        protein_readings[protein], protein.accession_number
                    )
                    residuals_all, r_squared_all = self._calcResidualsR2All(
                        protein_readings[protein]
                    )

                    metrics = {
                        "standard_deviation": std_dev,
                        **{f"{k}_average": v for k, v in base_metrics.items()},
                        "residuals_all": residuals_all,
                        "R_squared_all": r_squared_all,
                        "curve_fold_change": curve_fold_change,
                        "curve_peak": curve_peak,
                    }

                protein_metrics[protein] = metrics

            except Exception as e:
                print("Exception in _calculate_metrics")
                print(e)
                # TODO - take this out
                break

        return protein_metrics

    # TODO - this is straight up lifted from ICR, change and ideally use a library
    def _calcResidualsR2All(self, norm_abundances):
        """
        Calculate the residuals and the R squared for all the abundances from all the replicates for a protein.
        """

        # TODO - use generic stage names
        timepoint_map = {
            "Palbo": 0,
            "Late G1_1": 1,
            "G1/S": 2,
            "S": 3,
            "S/G2": 4,
            "G2_2": 5,
            "G2/M_1": 6,
            "M/Early G1": 7,
        }
        # if we have info for the protein in at least 2 replicates
        # TODO - repetitive
        # TODO - very similiar to _calcCurveFoldChange
        if len(norm_abundances) >= 2:
            reps = norm_abundances
            x = []
            for rep in reps:
                if rep != "Two":
                    continue
                for timepoint in norm_abundances[rep]:
                    x.append(timepoint_map[timepoint])
            x.sort()
            y = []
            for timepoint in timepoint_map:
                for rep in reps:
                    if rep != "Two":
                        continue
                    if timepoint in norm_abundances[rep]:
                        y.append(norm_abundances[rep][timepoint])

            p = np.poly1d(np.polyfit(x, y, 2))
            curve_abundances = p(x)
            residuals_all = np.polyfit(x, y, 2, full=True)[1][0]
            r_squared_all = round(r2_score(y, curve_abundances), 2)

        return residuals_all, r_squared_all

    # TODO - this is straight up lifted from ICR, change and ideally use a library
    def _calcCurveFoldChange(self, norm_abundances, uniprot_accession):
        """
        Calculates the curve_fold_change and curve peaks for the three or two replicates normalised abundance for each protein.
        """

        # TODO - use the stage_names for the project
        timepoint_map = {
            "Palbo": 0,
            "Late G1_1": 1,
            "G1/S": 2,
            "S": 3,
            "S/G2": 4,
            "G2_2": 5,
            "G2/M_1": 6,
            "M/Early G1": 7,
        }
        curves_all = {}

        # if we have info for the protein in at least 2 replicates
        # TODO - this code needs a tidy
        if len(norm_abundances) >= 2:
            reps = norm_abundances
            curves_all[uniprot_accession] = {}
            x = []
            for rep in reps:
                # TODO - make generic somehow
                if rep != "Two":
                    continue
                for timepoint in norm_abundances[rep]:
                    x.append(timepoint_map[timepoint])
            x.sort()
            y = []
            for timepoint in timepoint_map:
                for rep in reps:
                    # TODO - make generic somehow
                    if rep != "Two":
                        continue
                    if timepoint in norm_abundances[rep]:
                        y.append(norm_abundances[rep][timepoint])

            p = np.poly1d(np.polyfit(x, y, 2))
            curve_abundances = p(x)

            # find the timepoint peak of the curve
            curve_index = x[list(curve_abundances).index(max(curve_abundances))]
            for time_point, index in timepoint_map.items():
                if index == curve_index:
                    curve_peak = time_point

            # Calculate the fold change from the curve
            curve_fold_change = max(curve_abundances) / max(0.05, min(curve_abundances))

        return curve_fold_change, curve_peak

    # TODO - this is straight up lifted from ICR. Replace it, ideally with a library call
    # TODO - write a test for it first
    def _polyfit(self, x, y, degree):
        """
        Calculate the r-squared for polynomial curve fitting.
        """
        r_squared = 0
        coeffs = np.polyfit(x, y, degree)
        p = np.poly1d(coeffs)
        # calculate r-squared
        yhat = p(x)
        ybar = np.sum(y) / len(y)
        ssres = np.sum((y - yhat) ** 2)
        sstot = np.sum((y - ybar) ** 2)
        r_squared = ssres / sstot

        return 1 - round(r_squared, 2)

    def _impute(
        # TODO - all these pr types are wrong, and also probably bad variable names
        self,
        normalised_protein_readings: QuerySet[ProteinReading],
        replicates: QuerySet[Replicate],
        column_names: QuerySet[ColumnName],
    ):
        logger.info("impute missing data")

        replicates_by_name: dict = {}
        column_names_by_replicate: dict = {}

        for replicate in replicates:
            replicates_by_name[replicate.name] = replicate
            column_names_by_replicate[replicate.name] = []

        for column_name in column_names:
            column_names_by_replicate[column_name.replicate.name].append(
                column_name.sample_stage.name
            )

        imputed_protein_readings: dict = {}

        for protein in normalised_protein_readings.keys():
            imputed_protein_readings[protein] = {}

            for replicate_name in normalised_protein_readings[protein].keys():
                imputed_protein_readings[protein][replicate_name] = {}

                abundances_dict = normalised_protein_readings[protein][replicate_name]
                abundances = list(abundances_dict.values())

                stage_names = column_names_by_replicate[replicate_name]

                for idx, stage_name in enumerate(stage_names):
                    # Default value, should never be used
                    value = 0

                    if abundances_dict.get(stage_name, None) is not None:
                        value = abundances_dict[stage_name]
                    else:
                        last = None
                        next = None

                        # TODO - isn't there a better way to iterate?
                        for offset in range(1, len(stage_names)):
                            prev_idx = idx - offset
                            if prev_idx < 0:
                                # Gone before the beginning of the list, give up
                                break

                            prev_stage_name = stage_names[prev_idx]

                            if abundances_dict.get(prev_stage_name, None) is not None:
                                last = (offset, abundances_dict[prev_stage_name])
                                # last = abundances_dict[prev_stage_name]
                                break

                        for offset in range(1, len(abundances)):
                            # Look forward
                            # TODO - this seems to loop back to the beginning. Is that right?
                            next_idx = (idx + offset) % len(abundances)
                            next_stage_name = stage_names[next_idx]

                            if abundances_dict.get(stage_name, None) is not None:
                                next = (offset, abundances[next_stage_name])
                                # next = abundances[next_stage_name]
                                break

                        if last and next:
                            # Linear imputation between nearest timepoints
                            # TODO - find out why this calculation
                            # TODO - name variables better
                            d1, a1 = last
                            d2, a2 = next
                            step_height = (a1 - a2) / (d1 + d2)
                            value = d2 * step_height + a2

                    # imputed_protein_readings[protein][replicate_name][stage_name] = self._round(float(value))
                    # TODO - for some reason ICR rounds to 2, not 4. What to do?
                    imputed_protein_readings[protein][replicate_name][
                        stage_name
                    ] = round(float(value), 2)

        return imputed_protein_readings

    def _calculate_level_two_normalisation(
        self, normalised_protein_readings: QuerySet[ProteinReading], zero_min=False
    ):
        logger.info("relative log2 normalising abundances")

        level_two_normalised_protein_readings: dict = {}

        for protein in normalised_protein_readings:
            level_two_normalised_protein_readings[protein] = {}

            for replicate_name in normalised_protein_readings[protein]:
                level_two_normalised_protein_readings[protein][replicate_name] = {}

                min_value = 0
                # TODO - maybe this should be a function
                abundances = normalised_protein_readings[protein][replicate_name]

                abundance_values_non_null = [
                    val for val in abundances.values() if val is not None
                ]

                if not zero_min:
                    min_value = min(abundance_values_non_null)

                max_value = max(abundance_values_non_null)

                for stage_name, abundance in abundances.items():
                    denominator = max_value - min_value
                    if abundance is None or denominator == 0:
                        abundance_normalised = 0.5
                    else:
                        abundance_normalised = (abundance - min_value) / denominator

                    level_two_normalised_protein_readings[protein][replicate_name][
                        stage_name
                    ] = self._round(abundance_normalised)

        return level_two_normalised_protein_readings

    def _calculate_relative_log2_normalisation(
        self, normalised_protein_readings: QuerySet[ProteinReading], project: Project
    ):
        logger.info("relative log2 normalising abundances")

        log2_normalised_protein_readings: dict = {}

        log2_abundances: dict = {}

        for protein in normalised_protein_readings:
            log2_abundances[protein] = {}

            for replicate_name in normalised_protein_readings[protein]:
                log2_abundances[protein][replicate_name] = {}

                for stage_name in normalised_protein_readings[protein][replicate_name]:
                    log2_abundances[protein][replicate_name][stage_name] = math.log2(
                        normalised_protein_readings[protein][replicate_name][stage_name]
                    )

        for protein in log2_abundances:
            total_abundances = 0
            total_lengths = 0

            for replicate_name in log2_abundances[protein]:
                for stage_name in log2_abundances[protein][replicate_name]:
                    if log2_abundances[protein][replicate_name][stage_name] is not None:
                        total_abundances += log2_abundances[protein][replicate_name][
                            stage_name
                        ]
                        total_lengths += 1

            mean = total_abundances / total_lengths

            log2_normalised_protein_readings[protein] = {}
            for replicate_name in normalised_protein_readings[protein]:
                log2_normalised_protein_readings[protein][replicate_name] = {}

                for stage_name in normalised_protein_readings[protein][replicate_name]:
                    normalised_abundance = None

                    if log2_abundances[protein][replicate_name].get(stage_name):
                        normalised_abundance = self._round(
                            log2_abundances[protein][replicate_name][stage_name] - mean
                        )

                    log2_normalised_protein_readings[protein][replicate_name][
                        stage_name
                    ] = normalised_abundance

        return log2_normalised_protein_readings

    def _calculate_arrest_log2_normalisation(
        self, normalised_protein_readings: QuerySet[ProteinReading], project: Project
    ):
        logger.info("log2 arrest normalising abundances")

        log2_normalised_protein_readings: dict = {}

        # TODO - what should the stage name be?
        # TODO - is ARRESTING_AGENT the wrong name?
        ARRESTING_AGENT = "Nocodozole"

        # TODO - this is a hack, maybe add the field to the Project model?
        if project.name == "ICR":
            ARRESTING_AGENT = "Palbo"

        for protein in normalised_protein_readings.keys():
            log2_normalised_protein_readings[protein] = {}

            for replicate_name in normalised_protein_readings[protein]:
                if protein.accession_number == "Q09666":
                    print("+++ REPLICATES")

                    print(normalised_protein_readings[protein])
                if not log2_normalised_protein_readings[protein].get(replicate_name):
                    log2_normalised_protein_readings[protein][replicate_name] = {}

                for stage_name in normalised_protein_readings[protein][replicate_name]:
                    log2_reading = None
                    reading = normalised_protein_readings[protein][replicate_name][
                        stage_name
                    ]

                    if normalised_protein_readings[protein][replicate_name].get(
                        ARRESTING_AGENT
                    ):
                        arrest_reading = normalised_protein_readings[protein][
                            replicate_name
                        ][ARRESTING_AGENT]

                        if reading is not None and arrest_reading is not None:
                            log2_reading = self._round(
                                math.log2(reading / arrest_reading)
                            )

                    log2_normalised_protein_readings[protein][replicate_name][
                        stage_name
                    ] = log2_reading

        return log2_normalised_protein_readings

    def _calculate_first_level_normalisation(
        self, protein_readings: QuerySet[ProteinReading], medians
    ):
        logger.info("Normalising abundances")

        normalised_protein_readings: dict = {}

        protein_no = 0

        for protein_reading in protein_readings:
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Normalising for {protein_no}, {protein_reading.protein.accession_number} {protein_reading.column_name.sample_stage.name}",
            )

            reading = protein_reading.reading

            if reading:
                protein = protein_reading.protein
                replicate_name = protein_reading.column_name.replicate.name
                sample_stage_name = protein_reading.column_name.sample_stage.name

                median = medians[protein_reading.column_name.replicate][
                    protein_reading.column_name
                ]

                # TODO - need to round this?
                # TODO - might there be undefined medians? If so what to do?
                #   Maybe throw an exception in the median calculation func if one is encountered?
                #   They probably shouldn't happen.
                normalised_reading = reading / median

                if not normalised_protein_readings.get(protein):
                    normalised_protein_readings[protein] = {}

                if not normalised_protein_readings[protein].get(replicate_name):
                    normalised_protein_readings[protein][replicate_name] = {}

                normalised_protein_readings[protein][replicate_name][
                    sample_stage_name
                ] = normalised_reading

        return normalised_protein_readings

    def _calculate_means_across_replicates_by_stage(
        self,
        protein_readings: QuerySet[ProteinReading],
        with_bugs: bool,
        imputed: bool = False,
    ):
        logger.info("Calculating mean across replicates by stage for each protein")

        means: dict = {}

        protein_no = 0

        for protein in protein_readings.keys():
            protein_no += 1
            self._count_logger(
                protein_no,
                10000,
                f"Calculating mean for {protein_no}, {protein.accession_number}",
            )

            protein_reading = protein_readings[protein]
            abundances: dict = {}

            for replicate_name in protein_reading:
                for stage_name in protein_reading[replicate_name]:
                    if not abundances.get(stage_name):
                        abundances[stage_name] = []

                    if with_bugs and not imputed:
                        # We throw away the second reading
                        # TODO - how will this behave towards None?
                        if len(abundances[stage_name]) == 1:
                            continue

                    if protein_reading[replicate_name][stage_name] is not None:
                        abundances[stage_name].append(
                            protein_reading[replicate_name][stage_name]
                        )

            means[protein] = {}

            for stage_name in abundances:
                abundance = abundances[stage_name]

                if len(abundance):
                    mean = sum(abundance) / len(abundance)
                    means[protein][stage_name] = self._round(mean)
                else:
                    # TODO - is this the right thing to do?
                    means[protein][stage_name] = None

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
