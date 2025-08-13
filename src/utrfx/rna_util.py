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
        # self._fc_wt = self._fc_variant_generator()
    #     self._fc_variant = self._fc_variant_generator()

    # def _fc_wt_generator(self) -> fold_compound: 
    #     return fold_compound(self._fc_wt)
    
    # def _fc_variant_generator(self) -> fold_compound: 
    #     return fold_compound(self._fc_variant)
    
    def mfe_diff(self) -> float:
        fc_wt = RNA.fold_compound(self._wt_sequence)
        fc_variant = RNA.fold_compound(self._variant_sequence)

        wt_structure, wt_mfe = fc_wt.mfe()
        variant_structure, variant_mfe = fc_variant.mfe()

        return wt_mfe - variant_mfe
    
    def mfe_diff(self) -> float:            
        fc_wt = RNA.fold_compound(self._wt_sequence)
        fc_variant = RNA.fold_compound(self._variant_sequence)

        wt_structure, wt_mfe = fc_wt.mfe()
        variant_structure, variant_mfe = fc_variant.mfe()

        return wt_mfe - variant_mfe
