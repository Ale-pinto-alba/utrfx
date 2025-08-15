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
        self._pf()    

        wt_diversity = self._fc_wt.mean_bp_distance()
        variant_diversity = self._fc_variant.mean_bp_distance()

        return wt_diversity - variant_diversity
    
    def mfe_frequency_diff(self) -> float: 
        self._pf()    

        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()
        
        wt_mfe_frequency = self._fc_wt.pr_structure(wt_structure)
        variant_mfe_frequency = self._fc_variant.pr_structure(variant_structure)

        return wt_mfe_frequency - variant_mfe_frequency
    
    def _pf(self):
        self._fc_wt.pf()
        self._fc_variant.pf()