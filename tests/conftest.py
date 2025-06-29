import json
import pytest
from pathlib import Path
import pandas as pd

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
