import factory

from process.models import Project, Replicate


class ProjectFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = Project

    name = factory.Faker("word")
    proteome_file = factory.Faker("word")
    phosphoproteome_file = factory.Faker("word")


class ReplicateFactory(factory.django.DjangoModelFactory):
    class Meta:
        model = Replicate

    name = factory.Faker("word")
    project = factory.SubFactory(ProjectFactory)
    rank = factory.Sequence(lambda n: n + 1)
