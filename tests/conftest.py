import os
import pytest

from utrfx.genome import GenomeBuild, GRCh38

@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    return os.path.dirname(__file__)

@pytest.fixture(scope="session")
def fpath_data_dir(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data")

@pytest.fixture(scope="session")
def genome_build() -> GenomeBuild:
    return GRCh38