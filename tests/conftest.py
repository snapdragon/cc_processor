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

@pytest.fixture
def load_json():
    def _load(file_name):
        file_path = Path(__file__).parent / "fixtures" / file_name
        with open(file_path) as f:
            return json.load(f)
    return _load

@pytest.fixture
def load_pickle():
    def _load(file_name):
        file_path = Path(__file__).parent / "fixtures" / file_name

        return pd.read_pickle(file_path)

    return _load

@pytest.fixture
def basic_project_setup():
    project = Project.objects.get(name="SL")
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