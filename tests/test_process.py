import pytest
import numpy as np
from math import isclose

from process.factories import (
    ColumnNameFactory,
    ProjectFactory,
    ProteinFactory,
    ProteinReadingFactory,
    ReplicateFactory,
    SampleStageFactory,
)
from process.management.commands.process import Command
from process.models import ColumnName, ProteinReading


@pytest.mark.django_db
def test_all_replicates():
    command = Command()
    project = ProjectFactory()
    replicates = [
        ReplicateFactory(project=project, name="r1"),
        ReplicateFactory(project=project, name="r2"),
    ]

    def test_func(replicate_name):
        return {}

    results = command._all_replicates(func=test_func, replicates=replicates)

    assert results == {replicates[0].name: {}, replicates[1].name: {}}


@pytest.mark.django_db
def test_calculate_means():
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="r1")
    replicate2 = ReplicateFactory(name="r2", project=project)
    replicate3 = ReplicateFactory(name="r3", project=project)

    # TODO - why doesn't it care about these being different projects?
    #   More tests may be needed to check the code works with one of multiple projects.
    sample_stage1 = SampleStageFactory()
    sample_stage2 = SampleStageFactory()

    column_name1 = ColumnNameFactory(replicate=replicate1, sample_stage=sample_stage1)
    column_name2 = ColumnNameFactory(replicate=replicate2, sample_stage=sample_stage1)
    column_name3 = ColumnNameFactory(replicate=replicate3, sample_stage=sample_stage1)

    column_name4 = ColumnNameFactory(replicate=replicate1, sample_stage=sample_stage2)
    column_name5 = ColumnNameFactory(replicate=replicate2, sample_stage=sample_stage2)
    column_name6 = ColumnNameFactory(replicate=replicate3, sample_stage=sample_stage2)

    protein = ProteinFactory(project=project)
    protein_reading_1 = ProteinReadingFactory(
        protein=protein, column_name=column_name1, reading=1
    )
    protein_reading_2 = ProteinReadingFactory(
        protein=protein, column_name=column_name2, reading=4
    )
    protein_reading_3 = ProteinReadingFactory(
        protein=protein, column_name=column_name3, reading=7
    )
    protein_reading_4 = ProteinReadingFactory(
        protein=protein, column_name=column_name4, reading=1
    )
    protein_reading_5 = ProteinReadingFactory(
        protein=protein, column_name=column_name5, reading=None
    )
    protein_reading_6 = ProteinReadingFactory(
        protein=protein, column_name=column_name6, reading=5
    )

    readings = {
        replicate1.name: {
            sample_stage1.name: protein_reading_1.reading,
            sample_stage2.name: protein_reading_4.reading,
        },
        replicate2.name: {
            sample_stage1.name: protein_reading_2.reading,
            sample_stage2.name: protein_reading_5.reading,
        },
        replicate3.name: {
            sample_stage1.name: protein_reading_3.reading,
            sample_stage2.name: protein_reading_6.reading,
        },
    }

    results = command._calculate_means(readings, False)

    assert results == {
        sample_stage1.name: 4,
        sample_stage2.name: 3,
    }


@pytest.mark.django_db
def test_calc_replicate_column_medians():
    command = Command()

    project = ProjectFactory()

    replicate = ReplicateFactory(name="r1", project=project)

    # TODO - why doesn't it care about these being different projects?
    #   More tests may be needed to check the code works with one of multiple projects.
    sample_stage1 = SampleStageFactory()
    sample_stage2 = SampleStageFactory()

    column_name1 = ColumnNameFactory(replicate=replicate, sample_stage=sample_stage1)
    column_name2 = ColumnNameFactory(replicate=replicate, sample_stage=sample_stage2)

    protein = ProteinFactory(project=project)
    ProteinReadingFactory(protein=protein, column_name=column_name1, reading=1)
    ProteinReadingFactory(protein=protein, column_name=column_name1, reading=4)
    ProteinReadingFactory(protein=protein, column_name=column_name1, reading=4)
    ProteinReadingFactory(protein=protein, column_name=column_name2, reading=1)
    ProteinReadingFactory(protein=protein, column_name=column_name2, reading=None)
    ProteinReadingFactory(protein=protein, column_name=column_name2, reading=5)

    protein_readings = ProteinReading.objects.all()
    column_names = ColumnName.objects.all()

    results = command._calculate_replicate_stage_name_medians(
        replicate.name, protein_readings, column_names
    )

    assert results == {
        column_name1.sample_stage.name: 4.0,
        column_name2.sample_stage.name: 3.0,
    }


@pytest.mark.django_db
def test_tp():
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="rep_1")
    replicate2 = ReplicateFactory(name="rep_2", project=project)
    replicate3 = ReplicateFactory(name="rep_3", project=project)

    # TODO - why doesn't it care about these being different projects?
    #   More tests may be needed to check the code works with one of multiple projects.
    sample_stage1 = SampleStageFactory()
    sample_stage2 = SampleStageFactory()

    column_name1 = ColumnNameFactory(replicate=replicate1, sample_stage=sample_stage1)
    column_name2 = ColumnNameFactory(replicate=replicate2, sample_stage=sample_stage1)
    column_name3 = ColumnNameFactory(replicate=replicate3, sample_stage=sample_stage1)

    column_name4 = ColumnNameFactory(replicate=replicate1, sample_stage=sample_stage2)
    column_name5 = ColumnNameFactory(replicate=replicate2, sample_stage=sample_stage2)
    column_name6 = ColumnNameFactory(replicate=replicate3, sample_stage=sample_stage2)

    protein = ProteinFactory(project=project)
    protein_reading_1 = ProteinReadingFactory(
        protein=protein, column_name=column_name1, reading=1
    )
    protein_reading_2 = ProteinReadingFactory(
        protein=protein, column_name=column_name2, reading=4
    )
    protein_reading_3 = ProteinReadingFactory(
        protein=protein, column_name=column_name3, reading=7
    )
    protein_reading_4 = ProteinReadingFactory(
        protein=protein, column_name=column_name4, reading=1
    )
    protein_reading_5 = ProteinReadingFactory(
        protein=protein, column_name=column_name5, reading=None
    )
    protein_reading_6 = ProteinReadingFactory(
        protein=protein, column_name=column_name6, reading=5
    )

    readings = {
        replicate1.name: {
            sample_stage1.name: protein_reading_1.reading,
            sample_stage2.name: protein_reading_4.reading,
        },
        replicate2.name: {
            sample_stage1.name: protein_reading_2.reading,
            sample_stage2.name: protein_reading_5.reading,
        },
        replicate3.name: {
            sample_stage1.name: protein_reading_3.reading,
            sample_stage2.name: protein_reading_6.reading,
        },
    }

    sample_stage_1_result = command._tp(sample_stage1.name, readings)
    sample_stage_2_result = command._tp(sample_stage2.name, readings)

    assert sample_stage_1_result == [1.0, 4.0, 7.0]
    assert sample_stage_2_result == [1.0, 5.0]


@pytest.mark.django_db
def test_calculate_residuals_R2_all():
    command = Command()

    # TODO - this is just captured output, not reasoned. Fix.
    readings = {
        'One': {
            'Palbo': -0.2308, 'Late G1_1': 0.4059, 'G1/S': -0.7048, 'S': 0.0117, 'S/G2': -0.4986, 'G2_2': -0.3323, 'G2/M_1': -0.232, 'M/Early G1': 0.6562
        },
        'Two': {
            'Palbo': 0.1586, 'Late G1_1': 0.1537, 'G1/S': 0.1503, 'S': 0.0468, 'S/G2': 0.1054, 'G2_2': 0.1219, 'G2/M_1': 0.0836, 'M/Early G1': 0.1045
        }
    }

    result = command._calculate_residuals_R2_all(readings)

    assert result == (0.005632439523809521, 0.47)

@pytest.mark.django_db
def test_calculate_curve_fold_change():
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="One", project=project)
    replicate2 = ReplicateFactory(name="Two", project=project)

    readings = {
        'One': {
            'Palbo': -0.2308, 'Late G1_1': 0.4059, 'G1/S': -0.7048, 'S': 0.0117, 'S/G2': -0.4986, 'G2_2': -0.3323, 'G2/M_1': -0.232, 'M/Early G1': 0.6562
        },
        'Two': {
            'Palbo': 0.1586, 'Late G1_1': 0.1537, 'G1/S': 0.1503, 'S': 0.0468, 'S/G2': 0.1054, 'G2_2': 0.1219, 'G2/M_1': 0.0836, 'M/Early G1': 0.1045
        }
    }

    result = command._calculate_curve_fold_change(readings, [replicate1, replicate2])

    assert result == (1.8033637502236528, 'Palbo')


# TODO - write a deliberate, miniature version of this
# TODO - do one for SL as well
@pytest.mark.django_db
def test_calculate_protein_oscillation(load_json):
    command = Command()

    results_single = load_json("calculate_protein_oscillation_results_data.json")

    norm_method = "0-max"
    phosphosite = "212"

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="abundance_rep_1", project=project)
    replicate2 = ReplicateFactory(name="abundance_rep_2", project=project)

    result = command._calculate_protein_oscillation(
        results_single,
        norm_method,
        phosphosite,
        [replicate1, replicate2]
    )

    assert results_single["gene_name"] == "AHNAK"
    assert result == {
        "abundance_rep_1": {
            "Palbo": 0.45930000000000004,
            "Late G1_1": 0.14100000000000001,
            "G1/S": 0.5424,
            "S": 0.30889999999999995,
            "S/G2": 0.5092000000000001,
            "G2_2": 0.48460000000000003,
            "G2/M_1": 0.31079999999999997,
            "M/Early G1": -0.17710000000000004
        },
        "abundance_rep_2": {
            "Palbo": 0,
            "Late G1_1": -0.023700000000000054,
            "G1/S": -0.12629999999999997,
            "S": 0.027200000000000002,
            "S/G2": -0.03539999999999999,
            "G2_2": -0.03739999999999999,
            "G2/M_1": -0.10670000000000002,
            "M/Early G1": -0.1630999999999999
        }
    }

# TODO - find out why these results are what they are
# TODO - do one for SL as well
@pytest.mark.django_db
def test_calculate_metrics(load_json):
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="rep_1", project=project)
    replicate2 = ReplicateFactory(name="rep_2", project=project)

    sample_stage1 = SampleStageFactory(name='Palbo', rank=1, project=project)
    sample_stage2 = SampleStageFactory(name='Late G1_1', rank=2, project=project)
    sample_stage3 = SampleStageFactory(name='G1/S', rank=3, project=project)
    sample_stage4 = SampleStageFactory(name='S', rank=4, project=project)
    sample_stage5 = SampleStageFactory(name='S/G2', rank=5, project=project)
    sample_stage6 = SampleStageFactory(name='G2_2', rank=6, project=project)
    sample_stage7 = SampleStageFactory(name='G2/M_1', rank=7, project=project)
    sample_stage8 = SampleStageFactory(name='M/Early G1', rank=8, project=project)

    readings = {
        replicate1.name: {
            sample_stage1.name: -0.2308,
            sample_stage2.name: 0.4059,
            sample_stage3.name: -0.7048,
            sample_stage4.name: 0.0117,
            sample_stage5.name: -0.4986,
            sample_stage6.name: -0.3323,
            sample_stage7.name: -0.232,
            sample_stage8.name: 0.6562
        },
        replicate2.name: {
            sample_stage1.name: 0.1586,
            sample_stage2.name: 0.1537,
            sample_stage3.name: 0.1503,
            sample_stage4.name: 0.0468,
            sample_stage5.name: 0.1054,
            sample_stage6.name: 0.1219,
            sample_stage7.name: 0.0836,
            sample_stage8.name: 0.1045
        }
    }
        
    readings_average = {
        sample_stage1.name: -0.2308,
        sample_stage2.name: 0.4059,
        sample_stage3.name: -0.7048,
        sample_stage4.name: 0.0117,
        sample_stage5.name: -0.4986,
        sample_stage6.name: -0.3323,
        sample_stage7.name: -0.232,
        sample_stage8.name: 0.6562
    }

    metrics = command._calculate_metrics(
        readings,
        readings_average,
        [replicate1, replicate2],
        [sample_stage1, sample_stage2, sample_stage3, sample_stage4, sample_stage5, sample_stage6, sample_stage7, sample_stage8]
    )

    assert isclose(metrics['standard_deviation'], 0.3342441252911809, rel_tol=1e-5) == True
    assert isclose(metrics['variance_average'].item(), 0.18, rel_tol=1e-5) == True
    assert isclose(metrics['skewness_average'].item(), 0.04119630918711327, rel_tol=1e-5) == True
    assert isclose(metrics['kurtosis_average'].item(), 0.07170485201257046, rel_tol=1e-5) == True
    assert metrics['peak_average'] == 'M/Early G1'
    assert isclose(metrics['max_fold_change_average'], 1.361, rel_tol=1e-5) == True
    assert isclose(metrics['residuals_average'].item(), 0.9228990995833334, rel_tol=1e-5) == True
    assert isclose(metrics['R_squared_average'].item(), 0.36, rel_tol=1e-5) == True
    assert isclose(metrics['residuals_all'], 0.0056324395238095265, rel_tol=1e-5) == True
    assert isclose(metrics['R_squared_all'], 0.47, rel_tol=1e-5) == True
    assert isclose(metrics['curve_fold_change'], 1.8033637502236521, rel_tol=1e-5) == True
    assert metrics['curve_peak'] == 'Palbo'


# TODO - more of these
# TODO - do SL version
@pytest.mark.django_db
def test_calculate_ANOVA():
    command = Command()

    abundances = {
        'abundance_rep_1': {
            'Palbo': -0.2308,
            'Late G1_1': 0.4059,
            'G1/S': -0.7048,
            'S': 0.0117,
            'S/G2': -0.4986,
            'G2_2': -0.3323,
            'G2/M_1': -0.232,
            'M/Early G1': 0.6562
        },
        'abundance_rep_2': {
            'Palbo': 0.1586,
            'Late G1_1': 0.1537,
            'G1/S': 0.1503,
            'S': 0.0468,
            'S/G2': 0.1054,
            'G2_2': 0.1219,
            'G2/M_1': 0.0836,
            'M/Early G1': 0.1045
        }
    }

    anovas = command._calculate_ANOVA(abundances)

    assert anovas == {
        "p_value": 0.5788203095177293,
        "F_statistics": 0.849268808563235
    }


# TODO - more of these
# TODO - do SL version
@pytest.mark.django_db
def test_calculate_residuals_R2_all():
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="rep_1", project=project)
    replicate2 = ReplicateFactory(name="rep_2", project=project)

    readings = {
        'rep_1': {
            'Palbo': -0.2308,
            'Late G1_1': 0.4059,
            'G1/S': -0.7048,
            'S': 0.0117,
            'S/G2': -0.4986,
            'G2_2': -0.3323,
            'G2/M_1': -0.232,
            'M/Early G1': 0.6562
        },
        'rep_2': {
            'Palbo': 0.1586,
            'Late G1_1': 0.1537,
            'G1/S': 0.1503,
            'S': 0.0468,
            'S/G2': 0.1054,
            'G2_2': 0.1219,
            'G2/M_1': 0.0836,
            'M/Early G1': 0.1045
        }, # TODO - does the real version have 'abundance_average'?
        'abundance_average': {
            'Palbo': -0.2308,
            'Late G1_1': 0.4059,
            'G1/S': -0.7048,
            'S': 0.0117,
            'S/G2': -0.4986,
            'G2_2': -0.3323,
            'G2/M_1': -0.232,
            'M/Early G1': 0.6562
        }
    }

    residuals_all, r_squared_all = command._calculate_residuals_R2_all(readings, [replicate1, replicate2])


    assert isclose(residuals_all, 0.0056324395238095265, rel_tol=1e-5) == True
    assert isclose(r_squared_all, 0.47, rel_tol=1e-5) == True


# TODO - more of these
# TODO - do SL version
@pytest.mark.django_db
def test_generate_xs_ys():
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="rep_1", project=project)
    replicate2 = ReplicateFactory(name="rep_2", project=project)

    readings = {
        'rep_1': {
            'Palbo': -0.2308,
            'Late G1_1': 0.4059,
            'G1/S': -0.7048,
            'S': 0.0117,
            'S/G2': -0.4986,
            'G2_2': -0.3323,
            'G2/M_1': -0.232,
            'M/Early G1': 0.6562
        },
        'rep_2': {
            'Palbo': 0.1586,
            'Late G1_1': 0.1537,
            'G1/S': 0.1503,
            'S': 0.0468,
            'S/G2': 0.1054,
            'G2_2': 0.1219,
            'G2/M_1': 0.0836,
            'M/Early G1': 0.1045
        },
    }

    # TODO - add test for stage_names_map
    x, y, stage_names_map = command._generate_xs_ys(readings, [replicate1, replicate2])

    assert x == [0, 1, 2, 3, 4, 5, 6, 7]
    assert y == [0.1586, 0.1537, 0.1503, 0.0468, 0.1054, 0.1219, 0.0836, 0.1045]


# TODO - write a deliberate, miniature version of this
# TODO - do one for SL as well
@pytest.mark.django_db
def test_generate_phospho_regression_metrics(load_json):
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="One", project=project)
    replicate2 = ReplicateFactory(name="Two", project=project)

    sample_stage1 = SampleStageFactory(name='Palbo', rank=1, project=project)
    sample_stage2 = SampleStageFactory(name='Late G1_1', rank=2, project=project)
    sample_stage3 = SampleStageFactory(name='G1/S', rank=3, project=project)
    sample_stage4 = SampleStageFactory(name='S', rank=4, project=project)
    sample_stage5 = SampleStageFactory(name='S/G2', rank=5, project=project)
    sample_stage6 = SampleStageFactory(name='G2_2', rank=6, project=project)
    sample_stage7 = SampleStageFactory(name='G2/M_1', rank=7, project=project)
    sample_stage8 = SampleStageFactory(name='M/Early G1', rank=8, project=project)

    pre_phospho_metrics = load_json("pre_generate_phospho_metrics.json")
    post_phospho_metrics = load_json("post_generate_phospho_metrics.json")

    output = command._generate_phospho_regression_metrics(
        pre_phospho_metrics,
        [replicate1, replicate2],
        [sample_stage1, sample_stage2, sample_stage3, sample_stage4, sample_stage5, sample_stage6, sample_stage7, sample_stage8],
        with_bugs=True
    )

    assert output == post_phospho_metrics



# TODO - write a deliberate, miniature version of this
# TODO - do one for SL as well
@pytest.mark.django_db
def test_generate_protein_metrics(load_json):
    command = Command()

    project = ProjectFactory()

    replicate1 = ReplicateFactory(name="One", project=project)
    replicate2 = ReplicateFactory(name="Two", project=project)

    sample_stage1 = SampleStageFactory(name='Palbo', rank=1, project=project)
    sample_stage2 = SampleStageFactory(name='Late G1_1', rank=2, project=project)
    sample_stage3 = SampleStageFactory(name='G1/S', rank=3, project=project)
    sample_stage4 = SampleStageFactory(name='S', rank=4, project=project)
    sample_stage5 = SampleStageFactory(name='S/G2', rank=5, project=project)
    sample_stage6 = SampleStageFactory(name='G2_2', rank=6, project=project)
    sample_stage7 = SampleStageFactory(name='G2/M_1', rank=7, project=project)
    sample_stage8 = SampleStageFactory(name='M/Early G1', rank=8, project=project)

    pre_generate_protein_metrics = load_json("pre_generate_protein_metrics.json")
    post_generate_protein_metrics = load_json("post_generate_protein_metrics.json")

    output = command._generate_protein_metrics(
        pre_generate_protein_metrics,
        [replicate1, replicate2],
        [sample_stage1, sample_stage2, sample_stage3, sample_stage4, sample_stage5, sample_stage6, sample_stage7, sample_stage8],
        with_bugs=True
    )

    assert output == post_generate_protein_metrics


# TODO - write a deliberate, miniature version of this
# TODO - do one for SL as well
@pytest.mark.django_db
def test_generate_df(load_json):
    command = Command()

    pre_generate_df = load_json("pre_generate_df.json")
    post_generate_df = load_json("post_generate_df.json")

    output = command._generate_df(
        pre_generate_df,
    )

    assert output == post_generate_df
