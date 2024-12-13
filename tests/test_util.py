import pytest

from utrfx.model import FiveUTRCoordinates
from utrfx.util import download_fasta_from_ensembl, get_five_prime_sequence, uorf_extractor


@pytest.mark.online
def test_download_fasta_from_ensembl():
    fasta = download_fasta_from_ensembl("ENST00000381418")

    assert fasta.startswith("AGTT")

    assert len(fasta) > 0


@pytest.fixture(scope="module")
def five_utr_sequence(
    transcript_fasta: str, 
    five_utr: FiveUTRCoordinates,
) -> str:
    return get_five_prime_sequence(transcript_sequence=transcript_fasta, five_utrs=five_utr)


def test_uorf_extractor(
    five_utr: FiveUTRCoordinates, 
    five_utr_sequence: str,
):
    uorfs = uorf_extractor(five_utr=five_utr, five_sequence=five_utr_sequence)
    assert len(uorfs) == 3

    first_uorf, second_uorf, third_uorf = uorfs

    assert first_uorf.uorf.start == 16
    assert first_uorf.uorf.end == 67

    assert second_uorf.uorf.start == 302
    assert second_uorf.uorf.end == 407

    assert third_uorf.uorf.start == 510
    assert third_uorf.uorf.end == 576