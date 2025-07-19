import json
import pytest
from pathlib import Path
import pandas as pd

from process.factories import (
    ProjectFactory,
    ReplicateFactory,
    SampleStageFactory,
    ProteinFactory,
    PhosphoFactory,
)

from process.models import (
    Project,
    Replicate,
    SampleStage,
)

from process.constants import (
    PROJECT_SL
)

@pytest.fixture
def basic_project_setup():
    project = Project.objects.get(name=PROJECT_SL)
    replicates = Replicate.objects.filter(project=project, mean=False)
    sample_stages = SampleStage.objects.filter(project=project)

    # TODO - give proteins accession_numbers?
    proteins = [
        ProteinFactory(project=project),
        ProteinFactory(project=project),
        ProteinFactory(project=project),
    ]

    phosphos = [
        PhosphoFactory(protein=proteins[0]),
        PhosphoFactory(protein=proteins[0]),
        PhosphoFactory(protein=proteins[1]),
        PhosphoFactory(protein=proteins[1]),
        PhosphoFactory(protein=proteins[1]),
        PhosphoFactory(protein=proteins[2]),
        PhosphoFactory(protein=proteins[2]),
        PhosphoFactory(protein=proteins[2]),
        PhosphoFactory(protein=proteins[2]),
    ]

    return {
        "project": project,
        "replicates": replicates,
        "sample_stages": sample_stages,
        "proteins": proteins,
        "phosphos": phosphos,
    }

# TODO - get rid of this later
@pytest.fixture
def basic_project_setup_ICR():
    project = Project.objects.get(name="ICR")
    replicates = Replicate.objects.filter(project=project)
    sample_stages = SampleStage.objects.filter(project=project)

    # TODO - give proteins accession_numbers?
    proteins = [
        ProteinFactory(project=project),
        ProteinFactory(project=project),
        ProteinFactory(project=project),
    ]

    return {
        "project": project,
        "replicates": replicates,
        "sample_stages": sample_stages,
        "proteins": proteins,
    }