import pytest

from utrfx.model import FiveUTRCoordinates
from utrfx.util import fetch_fasta_from_ensembl, uorf_extractor


@pytest.mark.online
def test_fetch_fasta_from_ensembl_positive_strand():
     
     assert fetch_fasta_from_ensembl("ENST00000696628").startswith("GGTCGTTCCC") == True
     assert fetch_fasta_from_ensembl("ENST00000696628").endswith("CTATTTGAAA") == True


@pytest.mark.online
def test_fetch_fasta_from_ensembl_negative_strand():
     
     assert fetch_fasta_from_ensembl("ENST00000381418").startswith("AGTTGCGCTT") == True
     assert fetch_fasta_from_ensembl("ENST00000381418").endswith("ATAAGGGTAA") == True


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