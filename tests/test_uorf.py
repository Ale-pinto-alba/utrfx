import os

import pytest

from utrfx.genome import GenomeBuild, GenomicRegion, GRCh38, Strand
from utrfx.model import FiveUTR
from utrfx.uorf import download_fasta_from_ensembl, get_five_prime_sequence

@pytest.fixture(scope="module")
def transcript_fasta(fpath_test_dir: str) -> str:
    return download_fasta_from_ensembl("ENST00000381418", os.path.join(fpath_test_dir, "data", "tx_ENST00000381418"))

@pytest.fixture(scope="module")
def five_utr() -> FiveUTR:
    return FiveUTR([GenomicRegion(contig= GenomeBuild.contig_by_name(GRCh38, "8"), start= 22114419, end= 22115043, strand= Strand.NEGATIVE)])

@pytest.fixture(scope="module")
def five_utr_sequence(transcript_fasta: str, five_utr: FiveUTR) -> str:
    return get_five_prime_sequence(transcript_sequence= transcript_fasta, five_utrs= five_utr)

