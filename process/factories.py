import factory
from factory.django import DjangoModelFactory

from process.models import (
    ColumnName,
    Project,
    Protein,
    ProteinReading,
    Replicate,
    SampleStage,
)


class ProjectFactory(DjangoModelFactory):
    class Meta:
        model = Project

    name = factory.Faker("word")
    proteome_file = factory.Faker("word")
    phosphoproteome_file = factory.Faker("word")
    proteome_file_accession_number_column_name = factory.Faker("word")


class ReplicateFactory(DjangoModelFactory):
    class Meta:
        model = Replicate

    name = factory.Faker("word")
    project = factory.SubFactory(ProjectFactory)
    rank = factory.Sequence(lambda n: n + 1)


class SampleStageFactory(DjangoModelFactory):
    class Meta:
        model = SampleStage

    name = factory.Faker("word")
    rank = factory.Sequence(lambda n: n)
    project = factory.SubFactory(ProjectFactory)


class ColumnNameFactory(DjangoModelFactory):
    class Meta:
        model = ColumnName

    name = factory.Faker("word")
    sample_stage = factory.SubFactory(SampleStageFactory)
    replicate = factory.SubFactory(ReplicateFactory)


class ProteinFactory(DjangoModelFactory):
    class Meta:
        model = Protein

    accession_number = factory.Faker("word")


class ProteinReadingFactory(DjangoModelFactory):
    class Meta:
        model = ProteinReading

    column_name = factory.SubFactory(ColumnNameFactory)
    protein = factory.SubFactory(ProteinFactory)
    reading = factory.Faker("pyfloat", left_digits=2, right_digits=2, positive=True)
