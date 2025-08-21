import numpy as np

import RNA


class RNA_folding:
    """
    `RNA_folding` enables the extraction of features related to RNA folding and calculates the differences between
    a wild-type sequence and the one with a variant.

    For this class, it has been used the Python's API of the package ViennaRNA.
    See here: Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
    ViennaRNA Package 2.0
    Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
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
        self._wt_structure = self._wt_structure_generator()
        self._variant_structure = self._variant_structure_generator()

    def _fc_wt_generator(self) -> RNA.fold_compound:
        """"
        Return the `fold_compound` of the wild-type sequence. This object is the core data structure of the ViennaRNA package, that holds data related
        to the sequence, stores folding models and energy parameters, matrices and constraints.
        """
        return RNA.fold_compound(self._wt_sequence)
    
    def _fc_variant_generator(self) -> RNA.fold_compound: 
        """
        Return the `fold_compound`of the sequence with the variant.
        """
        return RNA.fold_compound(self._variant_sequence)
    
    def _wt_structure_generator(self) -> str:
        """
        Generate the structure as a dot-bracket representation of the wild-type sequence.
        """
        wt_structure, wt_mfe = self._fc_wt.mfe()
        return wt_structure
    
    def _variant_structure_generator(self) -> str:
        """
        Generate the structure as a dot-bracket representation of the sequence with the variant.
        """
        variant_structure, wt_mfe = self._fc_variant.mfe()
        return variant_structure
    
    def mfe_diff(self) -> float:
        """
        Calculate the difference between the Minimum Free Energy (MFE) ensemble of the wild-type and the altered sequence.

        Returns: an `float` with the difference.
        """
        wt_structure, wt_mfe = self._fc_wt.mfe()
        variant_structure, variant_mfe = self._fc_variant.mfe()
        return wt_mfe - variant_mfe
    
    def ensemble_diversity_diff(self) -> float:
        """
        Calculate the difference between the Ensemble diversity of the wild-type and the altered sequence.

        Returns: an `float` with the difference.
        """ 
        self._pf()    

        wt_diversity = self._fc_wt.mean_bp_distance()
        variant_diversity = self._fc_variant.mean_bp_distance()
        return wt_diversity - variant_diversity
    
    def mfe_frequency_diff(self) -> float: 
        """
        Give the difference between the MFE frequency of the wild-type and the altered sequence.

        Returns: an `float` with the difference.
        """ 
        self._pf()    

        wt_mfe_frequency = self._fc_wt.pr_structure(self._wt_structure)
        variant_mfe_frequency = self._fc_variant.pr_structure(self._variant_structure)
        return wt_mfe_frequency - variant_mfe_frequency
    
    def hamming_distance(self) -> int:         
        """
        Give the Hamming distance between the wild-type and the altered sequence. This can only be calculated if both sequences
        are the same length. 

        Returns: an `int` with the distance.
        """ 
        return RNA.hamming_distance(self._wt_structure, self._variant_structure)
    
    def bp_distance(self) -> int: 
        """
        Give the base pair distance between the wild-type and the altered sequence.

        Returns: an `int` with the distance.
        """ 
        return RNA.bp_distance(self._wt_structure, self._variant_structure)
    
    def _pf(self):
        """
        Calculate the partition function over all possible secondary structures for both sequences.
        """
        self._fc_wt.pf()
        self._fc_variant.pf()

    def unpaired_bases_diff(self) -> int: 
        """
        Calculate the difference between the number of unpaired bases between the two sequences.

        Returns: an `int` with the difference.
        """
        wt_unpaired = self._wt_structure.count('.')
        variant_unpaired = self._variant_structure.count('.')
        return abs(wt_unpaired - variant_unpaired)
    
    def unpaired_bases_percentage(self) -> float:
        """
        Calculate the percentage of unpaired bases in the altered sequence.

        Returns: a `float` with the percentage.
        """ 
        variant_unpaired = self._variant_structure.count('.')
        return (variant_unpaired * 100) / len(self._variant_structure)
    
    def number_loops_diff(self) -> int: 
        """
        Calculate the difference in the number of loops in the sequence with the variant.

        Returns: a `int` with the difference. 
        """  
        wt_pt = RNA.ptable(self._wt_structure)
        wt_loops = RNA.loopidx_from_ptable(wt_pt)
        variant_pt = RNA.ptable(self._variant_structure)
        variant_loops = RNA.loopidx_from_ptable(variant_pt)        

        wt_loop_idx = [wt_loops[i] for i in range(1, len(wt_pt))]
        wt_unique_loops = set(wt_loop_idx)
        var_loop_idx = [variant_loops[i] for i in range(1, len(variant_pt))]
        var_unique_loops = set(var_loop_idx)

        return len(wt_unique_loops) - len(var_unique_loops)
    
    def number_loops(self) -> int:
        """
        Count the number of loops in the sequence with the variant.

        Returns: a `int` with the number of loops. 
        """  
        variant_pt = RNA.ptable(self._variant_structure)
        variant_loops = RNA.loopidx_from_ptable(variant_pt)
        
        loop_idx = [variant_loops[i] for i in range(1, len(variant_pt))]
        unique_loops = set(loop_idx)

        return len(unique_loops)
    
    def structural_difference_at_variant_position(
        self,
        variant_pos: int,
    ) -> int:
        """
        Calculates the difference in structural depth between the sequences.

        Arg:
            variant_pos: `int` with the variant position within the sequence.

        Returns: an `int` indicating the difference in cumulative structure score at the variant position. 
        """ 
        def _dot_bracket_depth(structure: str, variant_pos: int) -> int:
            score = 0
            for i in range(variant_pos + 1):
                if structure[i] == '(':
                    score += 1
                elif structure[i] == ')':
                    score -= 1
            return score
        
        wt_score = _dot_bracket_depth(self._wt_structure, variant_pos)
        variant_score = _dot_bracket_depth(self._variant_structure, variant_pos)
        return wt_score - variant_score

    def variant_pos_structural_element(
        self,
        variant_pos: int,
    ) -> str: 
        """
        Indicate the type of structural element in the variant position of the altered sequence.

        Arg:
            variant_pos: `int` with the variant position within the sequence.

        Returns: a `str` with the element type.
        """
        if self._wt_structure[variant_pos] == "(" or self._wt_structure[variant_pos] == ")":
            return "Stem"
        else:
            return "Unpaired"

    @staticmethod
    def generate_ubox_probs(
        sequence: str, 
        threshold: float = 1e-5,
    ) -> list:
        """
        Stores the ubox base pairing probabilities between nucleotides of the sequence.

        Args:
            sequence: a `str` containing the nucleotide sequence.
            threshold: a `float` to filter base pair probabilities.

        Returns: a `list` containing the ubox base pair probabilities.
        """
        ubox = []
        fc = RNA.fold_compound(sequence)
        (propensity, ensemble_energy) = fc.pf()
        basepair_probs = fc.bpp()
        for i in range(1, len(sequence)+1):
            for j in range(i+1, len(sequence)+1):
                p = basepair_probs[i][j]
                if p > threshold:
                    ubox.append({"i": i, "j": j, "score": p})
        return ubox
    
    @staticmethod
    def generate_lbox_pairs(
        sequence: str, 
    ) -> list:
        """
        Stores the lbox pairs of the sequence.

        Args:
            sequence: a `str` containing the nucleotide sequence.

        Returns: a `list` containing the lbox pairs.
        """
        lbox = []
        fc = RNA.fold_compound(sequence)
        mfe_struct, mfe_energy = fc.mfe()
        stack = []
        mfe_pairs = set()
        for i, c in enumerate(mfe_struct):
            pos = i + 1
            if c == '(':
                stack.append(pos)
            elif c == ')':
                if stack:
                    j = stack.pop()
                mfe_pairs.add((min(j, pos), max(j, pos)))

        for pos1, pos2 in mfe_pairs:
            lbox.append({"pos1": pos1, "pos2": pos2, "score": 0.95})
        return lbox
    
    def compare_number_lbox_pairs(
        self,
        wt_lbox: list, 
        variant_lbox: list, 
    ) -> int:
        """
        Calculates the difference between the number of lbox pairs between both sequences.

        Args:
            lbox*: a `list` containing the lbox probabilities of a sequence.

        Returns: a `float` with the difference.
        """
        return len(wt_lbox) - len(variant_lbox)

    def compare_ubox_probs(
        self,
        wt_ubox: list, 
        variant_ubox: list, 
    ) -> float:
        """
        Compares and calculates the difference between all the ubox probabilities between both sequences.

        Args:
            ubox*: a `list` containing the ubox probabilities of a sequence.

        Returns: a `float` with the difference.
        """
        wt_sum_ubox = sum([x["score"] for x in wt_ubox])
        variant_sum_ubox = sum([x["score"] for x in variant_ubox])
        return wt_sum_ubox - variant_sum_ubox
    
    # def compare_ubox_pairs(
    #     self,
    #     wt_ubox: list, 
    #     variant_ubox: list, 
    # ) -> float:
    #     """
    #     Compares and calculates the difference between all the ubox probabilities between both sequences.

    #     Args:
    #         ubox*: a `list` containing the ubox probabilities of a sequence.

    #     Returns: a `float` with the difference.
    #     """
    #     return len(wt_ubox) - len(variant_ubox)