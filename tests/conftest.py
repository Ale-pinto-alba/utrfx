import os
import pytest

from utrfx.genome import GenomeBuild, GRCh38, GenomicRegion, Strand
from utrfx.uorf import UORFs

@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    return os.path.dirname(__file__)

@pytest.fixture(scope="session")
def fpath_data_dir(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data")

@pytest.fixture(scope="session")
def fpath_example_gtf(fpath_data_dir: str) -> str:
    return os.path.join(fpath_data_dir,  "Homo.sapiens.GRCh38_sample.gtf")

@pytest.fixture(scope="session")
def genome_build() -> GenomeBuild:
    return GRCh38

@pytest.fixture(scope="session")
def first_region_favorite_tx(genome_build: genome_build):
    contig = genome_build.contig_by_name("chr22")
    return GenomicRegion(contig=contig, start=44668712, end=44668805, strand=Strand.POSITIVE)

@pytest.fixture(scope="session")
def second_region_favorite_tx(genome_build: genome_build):
    contig = genome_build.contig_by_name("chr22")
    return GenomicRegion(contig=contig, start=44702491, end=44702501, strand=Strand.POSITIVE)

@pytest.fixture(scope="session")
def first_region_another_tx(genome_build:genome_build):
    contig = genome_build.contig_by_name("chr22")
    return GenomicRegion(contig=contig, start=26837999, end=26838057, strand=Strand.NEGATIVE)

@pytest.fixture(scope="session")
def second_region_another_tx(genome_build:genome_build):
    contig = genome_build.contig_by_name("chr22")
    return GenomicRegion(contig=contig, start=26841401, end=26841576, strand=Strand.NEGATIVE)

@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    return os.path.dirname(__file__)

@pytest.fixture(scope="session")
def example_uorfs(fpath_fasta: str):
    return UORFs(tx_id="ENST00000381418.9", five_prime_sequence=fpath_fasta)