import pytest

from utrfx.model import FiveUTRCoordinates
from utrfx.util import fetch_genomic_sequence_from_ensembl, uorf_extractor, get_five_prime_sequence

@pytest.mark.online
@pytest.mark.parametrize(
     "tx_id, start, end, n_bases",
     [
          ("ENST00000650287", "GTCAGACTCG", "TATTACCACA", 28_007), # tx on the + strand
          ("ENST00000381418", "AGTTGCGCTT", "ATAAGGGTAA", 16_592), # tx on the - strand
     ]
)
def test_fetch_genomic_sequence_from_ensembl(
     tx_id: str,
     start: str,
     end: str,
     n_bases: int,
):
     cdna = fetch_genomic_sequence_from_ensembl(tx_id)
     assert cdna.startswith(start)
     assert cdna.endswith(end)
     assert len(cdna) == n_bases


@pytest.mark.online
def test_get_five_prime_sequence(
     lpl_five_utr: FiveUTRCoordinates,
):
     genomic_sequence = fetch_genomic_sequence_from_ensembl("ENST00000650287") 
     five_utr_sequence = get_five_prime_sequence(genomic_sequence, lpl_five_utr) # tx on the + strand
     assert five_utr_sequence.startswith("GTCAGACTCG")
     assert five_utr_sequence.endswith("GCGCCCCGAG")
     assert len(five_utr_sequence) == 188


@pytest.mark.online
def test_get_five_prime_sequence(
     hr_five_utr: FiveUTRCoordinates,
):
     genomic_sequence = fetch_genomic_sequence_from_ensembl("ENST00000381418")
     five_utr_sequence = get_five_prime_sequence(genomic_sequence, hr_five_utr) # tx on the - strand
     assert five_utr_sequence.startswith("AGTTGCGCTT")
     assert five_utr_sequence.endswith("CAGGAGAGTG")
     assert len(five_utr_sequence) == 623


def test_uorf_extractor(
    hr_five_utr: FiveUTRCoordinates, 
    hr_five_utr_sequence: str,
):
    uorfs = uorf_extractor(five_utr=hr_five_utr, five_sequence=hr_five_utr_sequence)
    assert len(uorfs) == 3

    first_uorf, second_uorf, third_uorf = uorfs

    assert first_uorf.uorf.start == 16
    assert first_uorf.uorf.end == 67

    assert second_uorf.uorf.start == 302
    assert second_uorf.uorf.end == 407

    assert third_uorf.uorf.start == 510
    assert third_uorf.uorf.end == 576