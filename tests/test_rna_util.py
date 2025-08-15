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

    def test_ensemble_diversity_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.ensemble_diversity_diff()
        
        assert actual == pytest.approx(3, abs=0.1)
    
    def test_mfe_frequency_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.mfe_frequency_diff()
        
        assert actual == pytest.approx(-0.06, abs=0.1)

    def test_hamming_distance(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.hamming_distance()
        
        assert actual == 18

    def test_bp_distance(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.bp_distance()
        
        assert actual == 13
    
    def test_unpaired_bases_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.unpaired_bases_diff()
        
        assert actual == 10

    def test_variant_change(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.variant_pos_change_type(15)
        
        assert actual == False

    # def test_total_probs_diff(
    #     self,
    #     wt_sequence: str, 
    #     variant_sequence: str,
    # ):
    #     folding = RNA_folding(wt_sequence, variant_sequence)

    #     lbox_wt, ubox_wt = folding.generate_probs(wt_sequence)
    #     lbox_variant, ubox_variant = folding.generate_probs(variant_sequence)

    #     actual = folding.compare_probs(lbox_wt, ubox_wt, lbox_variant, ubox_variant)
        
    #     assert actual == 13