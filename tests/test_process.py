import pytest

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
def test_calculate_means_across_replicates_by_stage():
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

    results = command._calculate_means_across_replicates_by_stage(readings, False)

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
def test_calcResidualsR2All():
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

    result = command._calcResidualsR2All(readings)

    assert result == (0.005632439523809521, 0.47)