import factory
from factory.django import DjangoModelFactory

from process.models import (
    Project,
    Protein,
    Replicate,
    SampleStage,
    StatisticType,
    Statistic,
    Abundance,
    Phospho
)

class ProjectFactory(DjangoModelFactory):
    class Meta:
        model = Project

    name = factory.Faker("word")
    proteome_file = factory.Faker("word")
    phosphoproteome_file = factory.Faker("word")
    proteome_file_accession_number_column_name = factory.Faker("word")
    processable = False


class ReplicateFactory(DjangoModelFactory):
    class Meta:
        model = Replicate

    name = factory.Faker("word")
    project = factory.SubFactory(ProjectFactory)
    rank = factory.Sequence(lambda n: n + 1)
    mean = False


class SampleStageFactory(DjangoModelFactory):
    class Meta:
        model = SampleStage

    name = factory.Faker("word")
    rank = factory.Sequence(lambda n: n)
    project = factory.SubFactory(ProjectFactory)


class ProteinFactory(DjangoModelFactory):
    class Meta:
        model = Protein

    accession_number = factory.Faker("word")


class PhosphoFactory(DjangoModelFactory):
    class Meta:
        model = Phospho

    protein = factory.SubFactory(ProteinFactory)
    mod = factory.Faker("word")
    phosphosite = factory.Faker("word")



class StatisticTypeFactory(DjangoModelFactory):
    class Meta:
        model = StatisticType

    name = factory.Faker("word")


class StatisticFactory(DjangoModelFactory):
    class Meta:
        model = Statistic

    statistic_type = factory.SubFactory(StatisticTypeFactory)
    # Only one of these should be populated so no defaults are given here
    project = None
    protein = None
    phospho = None
    metrics = {}

class AbundanceFactory(DjangoModelFactory):
    class Meta:
        model = Abundance

    statistic = factory.SubFactory(StatisticFactory)
    replicate = factory.SubFactory(ReplicateFactory)
    sample_stage = factory.SubFactory(SampleStageFactory)
    reading = 1
