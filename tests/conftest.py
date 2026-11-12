import json
from pathlib import Path
import pytest

@pytest.fixture
def base_config():
    c=json.loads((Path(__file__).resolve().parents[1]/"examples/steady_plate.json").read_text())
    c["geometry"]["per_side"]=5
    c["geometry"]["radial"]=4
    return c
