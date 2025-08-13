from ViennaRNA import fold_compound

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
        self._fc_wt = self._fc_variant_generator()
        self._fc_variant = self._fc_variant_generator()

    def _fc_wt_generator(self) -> fold_compound: 
        return fold_compound(self._fc_wt)
    
    def _fc_variant_generator(self) -> fold_compound: 
        return fold_compound(self._fc_variant)
    
    def mfe_diff(self) -> int:
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()
        return (wt_mfe - variant_mfe)