import tempfile
import os
import subprocess

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
    
    def hamming_distance(self) -> float: 
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()
        
        return RNA.hamming_distance(wt_structure, variant_structure)
    
    def bp_distance(self) -> float: 
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()
        
        return RNA.bp_distance(wt_structure, variant_structure)
    
    def _pf(self):
        self._fc_wt.pf()
        self._fc_variant.pf()

    def unpaired_bases_diff(self) -> float: 
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()
        
        wt_unpaired = wt_structure.count('.')
        variant_unpaired = variant_structure.count('.')

        return abs(wt_unpaired - variant_unpaired)
    
    def variant_pos_change_structural_element(
        self,
        variant_pos: int,
    ) -> bool: 
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()

        if wt_structure[variant_pos] == variant_structure[variant_pos]:
            return True
        else:
            return False
        
    def variant_pos_structural_element(
        self,
        variant_pos: int,
    ) -> bool: 
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()

        if wt_structure[variant_pos] == "(" or wt_structure[variant_pos] == ")":
            return "Stem"
        else:
            return "Unpaired"
    

    # def compare_probs(
    #     self,
    #     lbox1: list, 
    #     ubox1: list, 
    #     lbox2: list, 
    #     ubox2: list,
    # ):
    #     sum_probs_first_seq = sum([x["score"] for x in lbox1 + ubox1])
    #     suma_probs_second_seq = sum([x["score"] for x in lbox2 + ubox2])
    #     probs_diff = sum_probs_first_seq - suma_probs_second_seq

    #     return probs_diff
    
    # def generate_probs(
    #     self,
    #     sequence: str,
    #     threshold: float = 1e-5
    # ):
    #     lbox = []
    #     ubox = []

    #     fc = RNA.fold_compound(sequence)
    #     fc.pf()  # Calcula la función de partición

    #     length = len(sequence)

    #     for i in range(1, length + 1):
    #         for j in range(i + 1, length + 1):
    #             prob = fc.pr_ij(i, j)  # 👈 corrección aquí
    #             if prob > threshold:
    #                 entry = {"pos1": i, "pos2": j, "score": prob}
    #                 if j - i == 1:
    #                     lbox.append(entry)
    #                 else:
    #                     ubox.append(entry)

    #     return lbox, ubox
    