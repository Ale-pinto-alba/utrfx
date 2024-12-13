import os
import pytest

from utrfx.genome import GenomicRegion, GenomeBuild, GRCh38, Strand
from utrfx.model import FiveUTRCoordinates


def pytest_addoption(parser):
    parser.addoption(
        "--runonline", action="store_true", default=False, help="run online tests"
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "online: mark test that require internet access to run"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runonline"):
        # --runonline given in cli: do not skip online tests
        return
    skip_online = pytest.mark.skip(reason="need --runonline option to run")
    for item in items:
        if "online" in item.keywords:
            item.add_marker(skip_online)


@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    return os.path.dirname(__file__)


@pytest.fixture(scope="session")
def fpath_data_dir(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data")


@pytest.fixture(scope="session")
def genome_build() -> GenomeBuild:
    return GRCh38


@pytest.fixture(scope="session")
def transcript_fasta() -> str:
    return "AGTTGCGCTTCTGGCGATGGCGATCAGAGGTCCTGCTGCGCTCTCCGCCGCGCTCTACCTCCATTAGCCGCGCTGCGCGGTGCTGCGCCCTCGCCGGTGCCTCTCTCCTGGGTCCCAGGATCGGCCCCCACCATCCAGGCACGACCCCCTTCCCCGGCCCCTCGGCCTTTCCCCCAACTCGGCCATCTCCGACCCGGGGCGCGTGTTCCCCCCGGCCCGGCGCCTTCTCTCCCTCCGGGGGCACCCGCTCCCTAGCCCCGGCCCGGCCCTCCCCGCGGCGCAGCACGGAGTCTCGGCGTCCCATGGCGCAACCTACGGCCTCGGCCCAGAAGCTGGTGCGGCCGATCCGCGCCGTGTGCCGCATCCTGCAGATCCCGGAGTCCGACCCCTCCAACCTGCGGCCCTAGAGCGCCCCCGCCGCCCCGGGGGAAGGAGAGCGCGAGCGCGCTGAGCAGACAGAGCGGGAGAACGCGTCCTCGCCCGCCGGCCGGGAGGCCCCGGAGCTGGCCCATGGGGAGCAGGCGCCCGGTGCCGGCCACGACGACCGCCACCGCCCGCGCCGCGACCGGCCGGTGAAGCCCAGGTAAGCGCCAGGAGCGCGCCGTCTGGGGACACTCGTGGCGG"


@pytest.fixture(scope="session")
def five_utr() -> FiveUTRCoordinates:
    return FiveUTRCoordinates([GenomicRegion(contig= GenomeBuild.contig_by_name(GRCh38, "8"), start= 22114419, end= 22115043, strand= Strand.NEGATIVE)])