import pytest
import numpy as np
from math import isclose
import pandas as pd
import pandas.testing as pdt
import math

from process.factories import (
    ProjectFactory,
    ProteinFactory,
    ReplicateFactory,
    SampleStageFactory,
    StatisticTypeFactory,
    StatisticFactory,
    AbundanceFactory,
)
from process.management.commands.process import Command
from process.models import Abundance, StatisticType, Statistic

from process.constants import (
    ABUNDANCES_RAW,
    PROTEIN_MEDIAN,
    ABUNDANCES_NORMALISED_MEDIAN,
    ABUNDANCES_NORMALISED_LOG2_ARREST,
    ABUNDANCES_NORMALISED_LOG2_MEAN,
    ABUNDANCES_NORMALISED_MIN_MAX,
    ABUNDANCES_NORMALISED_ZERO_MAX,
    PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN,
    PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN,
)

@pytest.mark.django_db
def test_calculate_medians(basic_project_setup):
    command = Command()

    project = basic_project_setup["project"]
    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    stat_type_prot_r = StatisticTypeFactory(name=ABUNDANCES_RAW)
    stat_1 = StatisticFactory(statistic_type=stat_type_prot_r, protein=proteins[0])
    stat_2 = StatisticFactory(statistic_type=stat_type_prot_r, protein=proteins[1])
    stat_3 = StatisticFactory(statistic_type=stat_type_prot_r, protein=proteins[2])

    reading = 0

    for stat in [stat_1, stat_2, stat_3]:
        for replicate in replicates:
            for sample_stage in sample_stages:
                reading += 1

                AbundanceFactory(statistic=stat, replicate=replicate, sample_stage=sample_stage, reading=reading)

    command._calculate_medians(
        project,
        replicates,
        sample_stages,
        True,
        with_bugs=False
    )

    stat_type_prot_median = StatisticType.objects.get(name=PROTEIN_MEDIAN)

    medians = Abundance.objects.filter(
        statistic__project=project,
        statistic__statistic_type=stat_type_prot_median
    )

    assert len(replicates) == 3
    assert len(sample_stages) == 10
    assert len(proteins) == 3
    assert len(medians) == 30

    median = 30

    for replicate in replicates:
        for sample_stage in sample_stages:
            median += 1

            abundance = Abundance.objects.get(statistic__project=project, replicate=replicate, sample_stage=sample_stage)

            assert abundance.reading == median


@pytest.mark.django_db
def test_calculate_normalised_medians(basic_project_setup):
    command = Command()

    project = basic_project_setup["project"]
    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    stat_type_prot_median, stat_prot_medians = create_readings(
        PROTEIN_MEDIAN,
        replicates,
        sample_stages,
        reading = 30,
        project=project
    )

    stat_type_prot_reading, stat_prot_readings = create_readings(
        ABUNDANCES_RAW,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    command._calculate_normalised_medians(proteins[0])        

    stat_type_normalised_median = StatisticType.objects.get(name=ABUNDANCES_NORMALISED_MEDIAN)
    stat_normalised_median = Statistic.objects.get(statistic_type=stat_type_normalised_median, protein=proteins[0])

    abundances = Abundance.objects.filter(statistic=stat_normalised_median)

    assert len(abundances) == 30

    reading = 0
    for abundance in abundances:
        reading += 1

        # Protein readings start at 1, protein medians at 31
        assert abundance.reading == reading / (reading + 30)



@pytest.mark.django_db
def test_calculate_means(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    _, stat_prot_readings = create_readings(
        ABUNDANCES_RAW,
        replicates,
        sample_stages,
        0,
        protein=proteins[0]
    )

    command._calculate_means(
        ABUNDANCES_RAW,
        proteins[0],
        False
    )

    abundances = Abundance.objects.filter(statistic=stat_prot_readings)

    assert len(abundances) == 40

    abundances = Abundance.objects.filter(
        statistic=stat_prot_readings,
        replicate__mean=True
    ).order_by(
        'sample_stage__rank'
    )

    assert len(abundances) == 10

    reading = 10

    # Readings are 1, 11, 21, then 12, 22, 32 and so on
    for abundance in abundances:
        reading += 1

        assert abundance.reading == reading


@pytest.mark.django_db
def test_calculate_arrest_log2_normalisation(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    stat_type_normalised_median, stat_normalised_median = create_readings(
        ABUNDANCES_NORMALISED_MEDIAN,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    command._calculate_arrest_log2_normalisation(
        proteins[0]
    )

    stat_type_arrest_log2_normalisation = StatisticType.objects.get(
        name=ABUNDANCES_NORMALISED_LOG2_ARREST
    )
    stat_arrest_log2_normalisation = Statistic.objects.get(
        statistic_type=stat_type_arrest_log2_normalisation,
        protein=proteins[0]
    )

    abundances = Abundance.objects.filter(statistic=stat_arrest_log2_normalisation)

    assert len(abundances) == 30

    assert len(abundances.filter(replicate__name="One")) == 10

    reading = 0

    # TODO - hardcoding "One" is bad
    for abundance in abundances.filter(
        replicate__name="One"
    ).order_by(
        'sample_stage__rank'
    ):
        reading += 1

        # The arresting agent is the last, hence divide by 10
        assert abundance.reading == round(math.log2(reading / 10), 4)


@pytest.mark.django_db
def test_calculate_relative_log2_normalisation(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    stat_type_normalised_median, stat_normalised_median = create_readings(
        ABUNDANCES_NORMALISED_MEDIAN,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    command._calculate_relative_log2_normalisation(proteins[0])

    stat_type_arrest_log2_mean = StatisticType.objects.get(
        name=ABUNDANCES_NORMALISED_LOG2_MEAN
    )
    stat_arrest_log2_mean = Statistic.objects.get(
        statistic_type=stat_type_arrest_log2_mean,
        protein=proteins[0]
    )

    abundances = Abundance.objects.filter(statistic=stat_arrest_log2_mean)

    print(abundances)

    sum = 0

    # 1 + 2 + 3 ... / 30
    for i in range(1, 31):
        sum += math.log2(i)

    mean = sum / 30

    assert len(abundances) == 30

    num = 0

    for abundance in abundances:
        num += 1

        assert abundance.reading == round(math.log2(num) - mean, 4)


@pytest.mark.django_db
def test_calculate_min_normalisation(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    stat_type_log2_mean, stat_log2_mean = create_readings(
        ABUNDANCES_NORMALISED_LOG2_MEAN,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    command._calculate_zero_or_min_normalisation(
        replicates,
        proteins[0],
    )

    stat_type_min_max = StatisticType.objects.get(
        name=ABUNDANCES_NORMALISED_MIN_MAX
    )
    stat_min_max = Statistic.objects.get(
        statistic_type=stat_type_min_max,
        protein=proteins[0]
    )

    abundances = Abundance.objects.filter(statistic=stat_min_max)

    # replicate one: min is 1, max is 10, so denominator is 9
    # replicate two: min is 11, max is 20, so denominator is 9
    for replicate in replicates:
        num = -1

        for ab in abundances.filter(replicate=replicate):
            num += 1

            assert round(ab.reading, 3) == round(num / 9, 3)

    assert len(abundances) == 30




@pytest.mark.django_db
def test_calculate_zero_max_normalisation(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    stat_type_normalised_median, stat_normalised_median = create_readings(
        ABUNDANCES_NORMALISED_MEDIAN,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    command._calculate_zero_or_min_normalisation(
        replicates,
        proteins[0],
        None,
        True,
    )

    stat_type_zero = StatisticType.objects.get(
        name=ABUNDANCES_NORMALISED_ZERO_MAX
    )
    stat_zero = Statistic.objects.get(
        statistic_type=stat_type_zero,
        protein=proteins[0]
    )

    return

    abundances = Abundance.objects.filter(statistic=stat_zero)

    # replicate one: min is 0, max is 10, so denominator is 10
    # replicate two: min is 0, max is 20, so denominator is 20
    replicate_num = 0
    num = 0

    for replicate in replicates:
        replicate_num += 1

        for ab in abundances.filter(replicate=replicate):
            num += 1

            assert round(ab.reading, 3) == round(num / (replicate_num * 10), 3)

    assert len(abundances) == 30



@pytest.mark.parametrize("statistic_type_name, phospho, phospho_ab, phospho_reg", [
    (ABUNDANCES_NORMALISED_LOG2_MEAN, False, False, False),
    (ABUNDANCES_NORMALISED_LOG2_MEAN, True, False, False),
    (PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN, True, True, False),
    (PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN, True, False, True),
])
@pytest.mark.django_db
def test_create_abundance_dataframe(statistic_type_name, phospho, phospho_ab, phospho_reg, basic_project_setup):
    command = Command()

    project = basic_project_setup["project"]
    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]
    phosphos = basic_project_setup["phosphos"]

    reading = 0

    if phospho:
        for ph in phosphos:
            create_readings(
                statistic_type_name,
                replicates,
                sample_stages,
                reading = reading,
                protein = None,
                phospho = ph
            )

            reading += 30
    else:
        for protein in proteins:
            create_readings(
                statistic_type_name,
                replicates,
                sample_stages,
                reading = reading,
                protein=protein
            )

            reading += 30

    df = command._create_abundance_dataframe(
        project,
        replicates,
        sample_stages,
        phospho,
        phospho_ab,
        phospho_reg
    )

    column_names = df.columns

    column_num = 0

    for sample_stage in sample_stages:
        for replicate in replicates:
            rep_stage_name = f"{replicate.name}_{sample_stage.name}"

            assert column_names[column_num] == rep_stage_name

            column_num += 1

    dec = 0
    count = 1
    extra = 0

    for index, row in df.iterrows():
        for col in df.columns:
            reading = (dec * 10) + count + extra

            assert float(row[col]) == reading

            if not reading % 30:
                extra += 20

            dec += 1

            if dec == 3:
                count += 1
                dec = 0


@pytest.mark.django_db
def test_calculate_metrics(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    statistic_type_name = ABUNDANCES_NORMALISED_LOG2_MEAN

    stat_type, stat = create_readings(
        statistic_type_name,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    command._calculate_means(
        statistic_type_name,
        protein = proteins[0],
        phospho = None,
        with_bugs = True
    )

    command._calculate_metrics(
        statistic_type_name,
        replicates,
        sample_stages,
        proteins[0],
    )

    stat = command._get_statistic(statistic_type_name, protein = proteins[0])

    metrics = stat.metrics

    numbers = {
        'R_squared_all': 1.0,
        'residuals_all': 8.313515829031342e-30,
        'kurtosis_average': 120.8625,
        'skewness_average': 0.0,
        'variance_average': 8.25,
        'R_squared_average': 1.0,
        'curve_fold_change': 1.8181818181818175, 
        'residuals_average': 2.5276769669211793e-30,
        'standard_deviation': 8.803408430829505,
        'max_fold_change_average': 9.0
    }

    for number in numbers:
        assert isclose(metrics[number], numbers[number]) == True

    assert metrics['curve_peak'] == 'Nocodozole'
    assert metrics['peak_average'] == 'Nocodozole'


@pytest.mark.django_db
def test_calculate_curve_fold_change(basic_project_setup):
    command = Command()

    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]

    _, stat = create_readings(
        ABUNDANCES_RAW,
        replicates,
        sample_stages,
        reading = 0,
        protein=proteins[0]
    )

    abundances = Abundance.objects.filter(statistic=stat)

    curve_fold_change, curve_peak = command._calculate_curve_fold_change(
        abundances,
        replicates,
        sample_stages,
    )

    print("++++++")
    print("++++++")
    print("++++++")
    print("++++++")
    print(curve_fold_change)
    print(curve_peak)

# def test_polyfit():
    



# TODO - unfinished
@pytest.mark.parametrize("statistic_type_name, phospho, phospho_ab, phospho_reg", [
    (ABUNDANCES_NORMALISED_LOG2_MEAN, False, False, False),
    # (ABUNDANCES_NORMALISED_LOG2_MEAN, True, False, False),
    # (PROTEIN_OSCILLATION_ABUNDANCES_LOG2_MEAN, True, True, False),
    # (PHOSPHO_REGRESSION_POSITION_ABUNDANCES_NORMALISED_LOG2_MEAN, True, False, True),
])
@pytest.mark.django_db
def test_calculate_fisher_g(statistic_type_name, phospho, phospho_ab, phospho_reg, basic_project_setup):
    command = Command()

    project = basic_project_setup["project"]
    replicates = basic_project_setup["replicates"]
    sample_stages = basic_project_setup["sample_stages"]
    proteins = basic_project_setup["proteins"]
    phosphos = basic_project_setup["phosphos"]

    reading = 0

    if phospho:
        for ph in phosphos:
            create_readings(
                statistic_type_name,
                replicates,
                sample_stages,
                reading = reading,
                protein = None,
                phospho = ph
            )

            reading += 30
    else:
        for protein in proteins:
            create_readings(
                statistic_type_name,
                replicates,
                sample_stages,
                reading = reading,
                protein=protein
            )

            reading += 30

    dict = command._calculate_fisher_g(
        project,
        replicates,
        sample_stages,
        phospho,
        phospho_ab,
        phospho_reg
    )

    # dict2 = command._calculate_fisher_g_r_version(
    #     project,
    #     replicates,
    #     sample_stages,
    #     phospho,
    #     phospho_ab,
    #     phospho_reg
    # )

    # print("+++++ DICT")
    # print("+++++ DICT")
    # print("+++++ DICT")
    # print("+++++ DICT")
    # print(dict)
    # print(dict2)



# @pytest.mark.django_db
# def test_generate_xs_ys(basic_project_setup):
#     command = Command()

#     replicates = basic_project_setup["replicates"]
#     sample_stages = basic_project_setup["sample_stages"]

#     readings = {}
#     reading = 50

#     for replicate in replicates:
#         reading -= 1

#         readings[replicate.name] = {}
#         for sample_stage in sample_stages:
#             readings[replicate.name][sample_stage.name] = reading

#     x, y, stage_names_map = command._generate_xs_ys(
#         replicates,
#         readings,
#         None,
#         True,
#     )

#     assert x == [0, 1, 2, 3, 4, 5, 6, 7]
#     assert y == [0.1586, 0.1537, 0.1503, 0.0468, 0.1054, 0.1219, 0.0836, 0.1045]




@pytest.mark.skip(reason="Broken")
@pytest.mark.django_db
def test_tp():
    command = Command()




# TODO - move this elsewhere
def create_readings(
    statistic_type_name, replicates, sample_stages, reading, project = None, protein = None, phospho=None
):
    stat_type = StatisticType.objects.get(name=statistic_type_name)

    if project:
        stat = StatisticFactory(statistic_type=stat_type, project=project)
    elif protein:
        stat = StatisticFactory(statistic_type=stat_type, protein=protein)
    elif phospho:
        stat = StatisticFactory(statistic_type=stat_type, phospho=phospho)
    else:
        raise Exception

    for replicate in replicates:
        for sample_stage in sample_stages:
            reading += 1

            Abundance.objects.create(statistic=stat, replicate=replicate, sample_stage=sample_stage, reading=reading)

    return stat_type, stat


