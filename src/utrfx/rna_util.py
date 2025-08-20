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
    
    def variant_pos_same_structural_element(
        self,
        variant_pos: int,
    ) -> bool:
        """
        Check if the type of structural element in the position of the variant is the same in both sequences, e.g. it goes
        from unpaired to stem.

        Returns: a `bool` indicating if there is a change (False) or not (True). 
        """ 
        if self._wt_structure[variant_pos] == self._variant_structure[variant_pos]:
            return True
        else:
            return False
        
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

    def generate_probs(
        self,    
        sequence: str, 
        threshold: float = 1e-5,
    ) -> list:
        """
        Stores the base pairing probabilities between nucleotides. First, the base pairing probability matrix is generated,
        and then sorted according to whether the elements are above the diagonal (ubox) or below the diagonal (lbox).

        Args:
            sequence: a `str` containing the nucleotide sequence.
            threshold: a `float` to filter base pair probabilities.

        Returns: two lists containing the base pair probabilities.
        """
        lbox = []
        ubox = []

        fc = RNA.fold_compound(sequence)
        fc.pf() 

        bpp_mat = fc.bpp()  

        for i in range(len(bpp_mat)):           
            row = bpp_mat[i]                    
            for j_offset, prob in enumerate(row):
                j = i + j_offset + 1             
                if prob > threshold:
                    pos1 = i + 1                
                    pos2 = j + 1
                    entry = {"pos1": pos1, "pos2": pos2, "score": prob}
                    if pos2 - pos1 == 1:
                        lbox.append(entry)
                    else:
                        ubox.append(entry)

        return lbox, ubox
    
    def compare_probs(
        self,
        lbox1: list, 
        ubox1: list, 
        lbox2: list, 
        ubox2: list,
    ) -> float:
        """
        Compares and calculates the difference between all the base pair probabilities between both sequences.

        Args:
            ubox*: a `list` containing the ubox probabilities of a sequence.
            lbox*: a `list` containing the lbox probabilities of a sequence.

        Returns: a `float` with the difference.
        """
        sum_probs_first_seq = sum([x["score"] for x in lbox1 + ubox1])
        suma_probs_second_seq = sum([x["score"] for x in lbox2 + ubox2])
        probs_diff = sum_probs_first_seq - suma_probs_second_seq

        return probs_diff
    