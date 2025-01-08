import pytest

from utrfx.model import FiveUTRCoordinates
from utrfx.util import download_fasta_from_ensembl, get_five_prime_sequence, uorf_extractor


@pytest.fixture(scope="module")
def five_utr_sequence(
    transcript_fasta: str, 
    hbb_five_utr: FiveUTRCoordinates,
) -> str:
    five_utr_sequence = get_five_prime_sequence(cdna_sequence=transcript_fasta, five_utrs=hbb_five_utr)

    assert len(five_utr_sequence) == 623

    return five_utr_sequence


@pytest.mark.online
def test_transcript_fasta(hbb_five_utr_sequence: str):

    assert download_fasta_from_ensembl("ENST00000381418").startswith("TTA") == True 

    assert download_fasta_from_ensembl("ENST00000381418").translate(str.maketrans("ATCG", "TAGC"))[::-1][:623] == hbb_five_utr_sequence
    # Check if the actual 5' UTR region (negative strand) matches the one downloaded from the ENSEMBL website (positive strand),
    # prior reverse complement


def test_uorf_extractor(
    hbb_five_utr: FiveUTRCoordinates, 
    hbb_five_utr_sequence: str,
):
    uorfs = uorf_extractor(five_utr=hbb_five_utr, five_sequence=hbb_five_utr_sequence)
    assert len(uorfs) == 3

    first_uorf, second_uorf, third_uorf = uorfs

    assert first_uorf.uorf.start == 16
    assert first_uorf.uorf.end == 67

    assert second_uorf.uorf.start == 302
    assert second_uorf.uorf.end == 407

    assert third_uorf.uorf.start == 510
    assert third_uorf.uorf.end == 576