import pytest

from process.factories import ProjectFactory, ReplicateFactory
from process.management.commands.process import Command


@pytest.mark.django_db
def test_all_replicates():
    project = ProjectFactory()
    replicates = [
        ReplicateFactory(project=project, name="r1"),
        ReplicateFactory(project=project, name="r2"),
    ]
    command = Command()

    def test_func(replicate):
        return {}

    results = command._all_replicates(func=test_func, replicates=replicates)

    assert results == {replicates[0]: {}, replicates[1]: {}}
