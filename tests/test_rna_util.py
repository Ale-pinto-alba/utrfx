import pytest

from utrfx.rna_util import RNA_folding

@pytest.fixture
def wt_sequence() -> str:
    return "ATACGATCATCGACGTACGCTACGATCATG"

@pytest.fixture
def variant_sequence() -> str:                           
    return "ATACGATCATCGACGAACGCTACGATCATG"

class TestRNAfolding:

    def test_variant_mfe(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.variant_mfe()
        
        assert actual == pytest.approx(-2, abs=0.1)

    def test_mfe_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.mfe_diff()
        
        assert actual == pytest.approx(-0.4, abs=0.1)

    def test_variant_ensemble_diversity(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.ensemble_diversity_diff()
        
        assert actual == pytest.approx(3, abs=0.1)

    def test_ensemble_diversity_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.ensemble_diversity_diff()
        
        assert actual == pytest.approx(3, abs=0.1)

    def test_variant_mfe_frequency(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.variant_mfe_frequency()
        
        assert actual == pytest.approx(0.2, abs=0.1)
    
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

    def test_unpaired_bases_percentage(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.unpaired_bases_percentage()
        
        assert actual == pytest.approx(73, abs=0.5)

    def test_number_loops_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.number_loops_diff()
        
        assert actual == 5

    def test_number_loops(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.number_loops()
        
        assert actual == 5

    @pytest.mark.parametrize(
        "variant_pos, expected",
        [
            (1, 0),    # Non changed position
            (15, 4),  # Position of the variant
        ]
    )
    def test_structural_difference_at_variant_position(
        self,
        wt_sequence: str, 
        variant_sequence: str,
        variant_pos: int,
        expected: int,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.structural_difference_at_variant_position(variant_pos)
        
        assert actual == expected

    @pytest.mark.parametrize(
        "variant_pos, expected",
        [
            (1, "Unpaired"),    
            (15, "Stem"),  
        ]
    )
    def test_variant_structural_element(
        self,
        wt_sequence: str, 
        variant_sequence: str,
        variant_pos: int,
        expected: bool,
    ):
        folding = RNA_folding(wt_sequence, variant_sequence)

        actual = folding.variant_pos_structural_element(variant_pos)
        
        assert actual == expected

    def test_lbox_probs_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        lbox_wt = RNA_folding.generate_lbox_pairs(wt_sequence)
        lbox_variant = RNA_folding.generate_lbox_pairs(variant_sequence)

        actual = RNA_folding.compare_number_lbox_pairs(lbox_wt, lbox_variant)

        assert actual == 5

    def test_jaccard_similarity(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        lbox_wt = RNA_folding.generate_lbox_pairs(wt_sequence)
        lbox_variant = RNA_folding.generate_lbox_pairs(variant_sequence)

        actual = RNA_folding.jaccard_similarity(lbox_wt, lbox_variant)

        assert actual == pytest.approx(0, abs=0.1)

    def test_ubox_total_probs_sum(
        self,
        variant_sequence: str,
    ):
        ubox_variant = RNA_folding.generate_ubox_probs(variant_sequence)

        actual = RNA_folding.ubox_total_probs_sum(ubox_variant)
        
        assert actual == pytest.approx(5.2, abs=0.1)

    def test_ubox_total_probs_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        ubox_wt = RNA_folding.generate_ubox_probs(wt_sequence)
        ubox_variant = RNA_folding.generate_ubox_probs(variant_sequence)

        actual = RNA_folding.ubox_total_probs_diff(ubox_wt, ubox_variant)
        
        assert actual == pytest.approx(0.99, abs=0.1)

    def test_ubox_mean_prob_diff(
        self,
        wt_sequence: str, 
        variant_sequence: str,
    ):
        ubox_wt = RNA_folding.generate_ubox_probs(wt_sequence)
        ubox_variant = RNA_folding.generate_ubox_probs(variant_sequence)

        actual = RNA_folding.ubox_mean_prob_diff(ubox_wt, ubox_variant)
        
        assert actual == pytest.approx(0.06, abs=0.1)

    def test_shannon_entropy(
        self,
        variant_sequence: str,
    ):
        ubox_variant = RNA_folding.generate_ubox_probs(variant_sequence)

        actual = RNA_folding.shannon_entropy(ubox_variant, len(variant_sequence))
        
        assert actual == pytest.approx(4.2, abs=0.1)

    def test_shannon_entropy_diff(
        self,
        wt_sequence: str,
        variant_sequence: str,
    ):
        ubox_wt = RNA_folding.generate_ubox_probs(wt_sequence)
        ubox_variant = RNA_folding.generate_ubox_probs(variant_sequence)

        wt_shannon_entropy = RNA_folding.shannon_entropy(ubox_wt, len(wt_sequence))        
        variant_shannon_entropy = RNA_folding.shannon_entropy(ubox_variant, len(variant_sequence))
        
        actual = RNA_folding.shannon_entropy_diff(wt_shannon_entropy, variant_shannon_entropy)

        assert actual == pytest.approx(0.4, abs=0.1)