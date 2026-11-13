"""Keep the package-version module and packaging metadata synchronized."""
import tomllib
from pathlib import Path

import thermofem2d


def test_internal_version_matches_package_metadata():
    pyproject=tomllib.loads((Path(__file__).resolve().parents[1]/'pyproject.toml').read_text())
    assert pyproject['project']['version'] == thermofem2d.__version__
