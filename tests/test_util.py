import pytest

from utrfx.model import FiveUTRCoordinates
from utrfx.util import fetch_genomic_sequence_from_ensembl, uorf_extractor


@pytest.mark.online
@pytest.mark.parametrize(
     "tx_id, start, end, n_bases",
     [
          ("ENST00000696628", "GGTCGTTCCC", "CTATTTGAAA", 2_412), # tx on the + strand
          ("ENST00000381418", "AGTTGCGCTT", "ATAAGGGTAA", 5_474), # tx on the - strand
     ]
)
def test_fetch_cdna_from_ensembl(
     tx_id: str,
     start: str,
     end: str,
     n_bases: int,
):
     cdna = fetch_genomic_sequence_from_ensembl(tx_id)
     assert cdna.startswith(start)
     assert cdna.endswith(end)
     assert len(cdna) == n_bases


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