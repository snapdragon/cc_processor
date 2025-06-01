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
from process.models import ProteinReading


@pytest.mark.django_db
def test_all_replicates():
    command = Command()
    project = ProjectFactory()
    replicates = [
        ReplicateFactory(project=project, name="r1"),
        ReplicateFactory(project=project, name="r2"),
    ]

    def test_func(replicate):
        return {}

    results = command._all_replicates(func=test_func, replicates=replicates)

    assert results == {replicates[0]: {}, replicates[1]: {}}


@pytest.mark.django_db
def test_calculate_abundance_means_across_replicates_by_stage():
    command = Command()

    replicate1 = ReplicateFactory(name="r1")
    replicate2 = ReplicateFactory(name="r1", project=replicate1.project)
    replicate3 = ReplicateFactory(name="r1", project=replicate1.project)

    sample_stage1 = SampleStageFactory()
    sample_stage2 = SampleStageFactory()

    column_name1 = ColumnNameFactory(replicate=replicate1, sample_stage=sample_stage1)
    column_name2 = ColumnNameFactory(replicate=replicate2, sample_stage=sample_stage1)
    column_name3 = ColumnNameFactory(replicate=replicate3, sample_stage=sample_stage1)

    column_name4 = ColumnNameFactory(replicate=replicate1, sample_stage=sample_stage2)
    column_name5 = ColumnNameFactory(replicate=replicate2, sample_stage=sample_stage2)
    column_name6 = ColumnNameFactory(replicate=replicate3, sample_stage=sample_stage2)

    protein = ProteinFactory()
    ProteinReadingFactory(protein=protein, column_name=column_name1, reading=1)
    ProteinReadingFactory(protein=protein, column_name=column_name2, reading=4)
    ProteinReadingFactory(protein=protein, column_name=column_name3, reading=7)
    ProteinReadingFactory(protein=protein, column_name=column_name4, reading=1)
    ProteinReadingFactory(protein=protein, column_name=column_name5, reading=None)
    ProteinReadingFactory(protein=protein, column_name=column_name6, reading=5)

    protein_readings = ProteinReading.objects.all()

    results = command._calculate_abundance_means_across_replicates_by_stage(
        protein_readings
    )

    assert results == {
        protein: {
            sample_stage1: 4,
            sample_stage2: 3,
        }
    }
