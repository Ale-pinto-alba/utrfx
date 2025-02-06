import os
import pytest

from utrfx.genome import GenomicRegion, GenomeBuild, GRCh38, Strand, VariantCoordinates
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
def hr_five_utr(
    genome_build: GenomeBuild
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic region corresponding to one of the transcripts of the HR gene (ENSEMBL transcript ID: `ENST00000381418.9`).

    Both Genomic Regions were obtained from the chromosome 8 GTF file.

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
def hr_five_utr_sequence() -> str:
    """
    5'UTR cDNA sequence of the transcript of the HR gene (ENSEMBL transcript ID: `ENST00000381418.9`) taken directly from
    the ENSEMBL website.
    
    see here: https://www.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000168453;r=8:22114419-22133384;t=ENST00000381418.
    """
    return "AGTTGCGCTTCTGGCGATGGCGATCAGAGGTCCTGCTGCGCTCTCCGCCG" \
        + "CGCTCTACCTCCATTAGCCGCGCTGCGCGGTGCTGCGCCCTCGCCGGTGC" \
        + "CTCTCTCCTGGGTCCCAGGATCGGCCCCCACCATCCAGGCACGACCCCCT" \
        + "TCCCCGGCCCCTCGGCCTTTCCCCCAACTCGGCCATCTCCGACCCGGGGC" \
        + "GCGTGTTCCCCCCGGCCCGGCGCCTTCTCTCCCTCCGGGGGCACCCGCTC" \
        + "CCTAGCCCCGGCCCGGCCCTCCCCGCGGCGCAGCACGGAGTCTCGGCGTC" \
        + "CCATGGCGCAACCTACGGCCTCGGCCCAGAAGCTGGTGCGGCCGATCCGC" \
        + "GCCGTGTGCCGCATCCTGCAGATCCCGGAGTCCGACCCCTCCAACCTGCG" \
        + "GCCCTAGAGCGCCCCCGCCGCCCCGGGGGAAGGAGAGCGCGAGCGCGCTG" \
        + "AGCAGACAGAGCGGGAGAACGCGTCCTCGCCCGCCGGCCGGGAGGCCCCG" \
        + "GAGCTGGCCCATGGGGAGCAGGCGCCCGGTGCCGGCCACGACGACCGCCA" \
        + "CCGCCCGCGCCGCGACCGGCCGGTGAAGCCCAGGGACCCCCCTCTGGGAG" \
        + "AGCCCCATGAGGGCAGGAGAGTG"


@pytest.fixture(scope="session")
def lpl_five_utr(
    genome_build = GenomeBuild,
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic region corresponding to one of the transcripts of the LPL gene (ENSEMBL transcript ID: `ENST00000650287.1`).

    Genomic Region was obtained from the chromosome 8 GTF file.
    """
    contig = genome_build.contig_by_name("8")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=19_939_252,
                end=19_939_440,
                strand=Strand.POSITIVE
                ),
        )
    ) 


@pytest.fixture(scope="session")
def hr_variant(
    genome_build: GenomeBuild,
) -> VariantCoordinates:
    """
    Single nucleotide variant (8-22130606-A-G) of the HR gene.
    """
    return VariantCoordinates.from_vcf_literal(contig=genome_build.contig_by_name("8"), pos= 22_130_606, ref= "A", alt= "G")


@pytest.fixture(scope="session")
def adsl_five_utr(
    genome_build: GenomeBuild,
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic region corresponding to one of the transcripts of the ADSL gene (ENSEMBL transcript ID: `ENST00000623063.3`).

    Genomic Region was obtained from the chromosome 22 GTF file.
    """
    contig = genome_build.contig_by_name("22")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=40_346_499,
                end=40_346_558,
                strand=Strand.POSITIVE
                ),
        )
    )


@pytest.fixture(scope="session")
def adsl_variant(
    genome_build: GenomeBuild,
) -> VariantCoordinates:
    """
    Single nucleotide variant (22-40346510-T-C) of the ADSL gene.
    """
    return VariantCoordinates.from_vcf_literal(contig=genome_build.contig_by_name("22"), pos= 40_346_510, ref= "T", alt= "C")


@pytest.fixture(scope="session")
def msh2_five_utr(
    genome_build: GenomeBuild,
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic region corresponding to one of the transcripts of the MSH2 gene (ENSEMBL transcript ID: `ENST00000543555.6`).

    Both Genomic Regions were obtained from the chromosome 2 GTF file.
    """
    contig = genome_build.contig_by_name("2")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=47_403_066,
                end=47_403_175,
                strand=Strand.POSITIVE
                ),
            GenomicRegion(
                contig=contig,
                start=47_403_359,
                end=47_403_389,
                strand=Strand.POSITIVE
                ),
        )
    )


@pytest.fixture(scope="session")
def msh2_variant(
    genome_build: GenomeBuild,
) -> VariantCoordinates:
    """
    Insertion variant (2-47403110-G-GA) of the MSH2 gene.
    """
    return VariantCoordinates.from_vcf_literal(contig=genome_build.contig_by_name("2"), pos= 47_403_110, ref= "G", alt= "GA")


@pytest.fixture(scope="session")
def lpl_five_utr(
    genome_build: GenomeBuild,
) -> FiveUTRCoordinates:
    """
    5'UTR Genomic region corresponding to one of the transcripts of the LPL gene (ENSEMBL transcript ID: `ENST00000650287.1`).

    Both Genomic Regions were obtained from the chromosome 8 GTF file.
    """
    contig = genome_build.contig_by_name("8")
    assert contig is not None

    return FiveUTRCoordinates(
        regions=(
            GenomicRegion(
                contig=contig,
                start=19_939_253,
                end=19_939_440,
                strand=Strand.POSITIVE
                ),
        )
    )


@pytest.fixture(scope="session")
def lpl_variant(
    genome_build: GenomeBuild,
) -> VariantCoordinates:
    """
    Single nucleotide variant (8-19939160-T-G) of the LPL gene.
    """
    return VariantCoordinates.from_vcf_literal(contig=genome_build.contig_by_name("8"), pos= 19_939_160, ref= "T", alt= "G")


@pytest.fixture(scope="session")
def lpl_fake_variant(
    genome_build: GenomeBuild,
) -> VariantCoordinates:
    """
    Fake single variant of the LPL gene.
    """
    return VariantCoordinates.from_vcf_literal(contig=genome_build.contig_by_name("8"), pos= 19_939_320, ref= "A", alt= "C")