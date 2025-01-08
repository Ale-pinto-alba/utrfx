import os
import pytest

from utrfx.genome import GenomicRegion, GenomeBuild, GRCh38, Strand
from utrfx.model import FiveUTRCoordinates
from utrfx.util import download_fasta_from_ensembl


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
def hbb_five_utr(
    genome_build: GenomeBuild
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic region corresponding to one of the transcripts of the HR gene (ENSEMBL transcript ID: `ENST00000381418.9`).

    Start coordinate was obtained directly from the ENSEMBL website and the end coordinate was calculated taking into account
    the 5'UTR length (cDNA sequence, therefore, number of bases, available on the website).

    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000168453;r=8:22114419-22133384;t=ENST00000381418
    """
    contig = genome_build.contig_by_name("8")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=22_130_427,
                end=22_131_010,
                strand=Strand.POSITIVE
                ).with_strand(other=Strand.NEGATIVE),
            GenomicRegion(
                contig=contig,
                start=22_129_170,
                end=22_129_210,
                strand=Strand.POSITIVE
                ).with_strand(other=Strand.NEGATIVE),
        )
    )


@pytest.fixture(scope="session")
@pytest.mark.online
def transcript_fasta():
    fasta = download_fasta_from_ensembl("ENST00000381418")

    assert fasta.startswith("AGTT")

    assert len(fasta) > 0

    return fasta


@pytest.fixture(scope="session")
def hbb_five_utr_sequence() -> str:
    """
    5'UTR cDNA sequence of the transcript of the HR gene (ENSEMBL transcript ID: `ENST00000381418.9`) taken directly from
    the ENSEMBL website.
    
    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000168453;r=8:22114419-22133384;t=ENST00000381418.
    """
    return "AGTTGCGCTTCTGGCGATGGCGATCAGAGGTCCTGCTGCGCTCTCCGCCGCGCTCTACCTCCATTAGCCGCGCTGCGCGGTGCTGCGCCCTCGCCGGTGCCTCTCTCCTGGGTCCCAGGATCGGCCCCCACCATCCAGGCACGACCCCCTTCCCCGGCCCCTCGGCCTTTCCCCCAACTCGGCCATCTCCGACCCGGGGCGCGTGTTCCCCCCGGCCCGGCGCCTTCTCTCCCTCCGGGGGCACCCGCTCCCTAGCCCCGGCCCGGCCCTCCCCGCGGCGCAGCACGGAGTCTCGGCGTCCCATGGCGCAACCTACGGCCTCGGCCCAGAAGCTGGTGCGGCCGATCCGCGCCGTGTGCCGCATCCTGCAGATCCCGGAGTCCGACCCCTCCAACCTGCGGCCCTAGAGCGCCCCCGCCGCCCCGGGGGAAGGAGAGCGCGAGCGCGCTGAGCAGACAGAGCGGGAGAACGCGTCCTCGCCCGCCGGCCGGGAGGCCCCGGAGCTGGCCCATGGGGAGCAGGCGCCCGGTGCCGGCCACGACGACCGCCACCGCCCGCGCCGCGACCGGCCGGTGAAGCCCAGGTAAGCGCCAGGAGCGCGCCGTCTGGGGACACTCGTGGCG"