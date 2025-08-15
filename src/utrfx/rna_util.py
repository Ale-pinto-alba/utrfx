import RNA

class RNA_folding:
    """
    """
    def __init__(
        self,
        wt_sequence: str,
        variant_sequence: str,
    ):
        self._wt_sequence = wt_sequence
        self._variant_sequence = variant_sequence
        self._fc_wt = self._fc_wt_generator()
        self._fc_variant = self._fc_variant_generator()

    def _fc_wt_generator(self) -> RNA.fold_compound: 
        return RNA.fold_compound(self._wt_sequence)
    
    def _fc_variant_generator(self) -> RNA.fold_compound: 
        return RNA.fold_compound(self._variant_sequence)
    
    def mfe_diff(self) -> float:
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()

        return wt_mfe - variant_mfe
    
    def ensemble_diversity_diff(self) -> float: 
        self._fc_wt.pf()
        self._fc_variant.pf()           

        wt_diversity = self._fc_wt.mean_bp_distance()
        variant_diversity = self._fc_variant.mean_bp_distance()

        return wt_diversity - variant_diversity