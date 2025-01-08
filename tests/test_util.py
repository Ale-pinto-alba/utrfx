import pytest

from utrfx.model import FiveUTRCoordinates
from utrfx.util import get_five_prime_sequence, uorf_extractor


@pytest.fixture(scope="module")
def five_utr_sequence(
    transcript_fasta: str, 
    hbb_five_utr: FiveUTRCoordinates,
) -> str:
    five_utr_sequence = get_five_prime_sequence(transcript_sequence=transcript_fasta, five_utrs=hbb_five_utr)

    assert len(five_utr_sequence) == 623

    return five_utr_sequence


def test_uorf_extractor(
    hbb_five_utr: FiveUTRCoordinates, 
    five_utr_sequence: str,
):
    uorfs = uorf_extractor(five_utr=hbb_five_utr, five_sequence=five_utr_sequence)
    assert len(uorfs) == 3

    first_uorf, second_uorf, third_uorf = uorfs

    assert first_uorf.uorf.start == 16
    assert first_uorf.uorf.end == 67

    assert second_uorf.uorf.start == 302
    assert second_uorf.uorf.end == 407

    assert third_uorf.uorf.start == 510
    assert third_uorf.uorf.end == 576