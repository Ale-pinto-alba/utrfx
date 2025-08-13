import pytest

from utrfx.rna_util import RNA_folding

@pytest.fixture
def wt_sequence() -> str:
    return "ATACGATCATCGACGTACGCTACGATCATG"

@pytest.fixture
def variant_sequence() -> str:                           
    return "ATACGATCATCGACGAACGCTACGATCATG"

class TestRNAfolding:

    def test_mfe_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.mfe_diff()
        
        assert actual == pytest.approx(-0.4, abs=0.1)