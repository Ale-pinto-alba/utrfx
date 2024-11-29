import typing

class UORFs:
    """
    `UORFs` represents the uORFs available in a 5'UTR region of a transcript (no uORFs whose end overlaps with the mORF).

    :param tx_id: GENCODE transcript identifier.
    :param five_prime_sequence: string or file containing the 5'UTR cDNA sequence.
    """
    def __init__(self, tx_id: str, five_prime_sequence: str):
        self._tx_id = tx_id
        self._five_prime_sequence = five_prime_sequence

    @property
    def tx_id(self) -> str:
        return self._tx_id

    @property
    def five_prime_sequence(self) -> str:
        if self._five_prime_sequence == str:
            return self._five_prime_sequence
        else:
            with open(self._five_prime_sequence) as f: 
                sequence = f.read().replace('\n', '')
            return sequence

    def uorf_extractor(self) -> typing.Collection[str]:
        """
        Take the nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.
        """
        five_sequence = self._five_prime_sequence
        uorfs = []
            
        while "ATG" in five_sequence:  
            codon_list = []  

            start_index = five_sequence.find("ATG")
            assert start_index != -1, "No ATG in the 5'UTR sequence."

            found_stop = False 
                
            for i in range(start_index, len(five_sequence) - len(five_sequence) % 3, 3):
                codon = five_sequence[i:i + 3]
                assert len(codon) == 3, "Codons must be 3nt long."

                codon_list.append(codon)

                if codon in ["TAA", "TAG", "TGA"]:
                    found_stop = True
                    break

            if found_stop:
                uorfs.append("".join(codon_list))

            five_sequence = five_sequence[start_index + 3:]

        return uorfs
    
    @property
    def uorfs(self) -> typing.Collection[str]:
        return self.uorf_extractor()
    
    def __repr__(self):
        return f"UORFS(tx_id= {self._tx_id}, uORFs{self.uorfs})"
    

class UORFs_calculations:
    """
    `UORF_calculations` is a container for the inner calculations of the uORFs.

    :param uorfs: collection of uORF(s) cDNA sequence(s) of a transcript.
    """
    def __init__(self, uorfs: UORFs):
        self._uorfs = uorfs
    
    def number_uorfs(self) -> int:
        return len(self._uorfs.uorfs)
    
    def uorfs_lengths(self) -> typing.Collection[int]:
        return [len(uorf) for uorf in self._uorfs.uorfs]
    
    def gc_content(self) -> typing.Collection[float]:
        """
        Get the GC content of each uORF.
        """
        gc_content_per_uorf = []

        for uorf in self._uorfs.uorfs:
            total = len(uorf)
            g = uorf.count("G")
            c = uorf.count("C")
            
            gc_content = ((g+c) / total) * 100
            gc_content_per_uorf.append(gc_content)

        return gc_content_per_uorf
    
    def __repr__(self):
        return f"UORFs_calculations(tx_id= {self._uorfs.tx_id}, uORFs= {self.number_uorfs()}, lengths= {self.uorfs_lengths()}, GC_content= {self.gc_content()})"


class UORFs_indel_analysis:
    """
    `UORFs_indel_analysis` is a container for the uORFs and the twenty nucleotides after (if possible) for their use
    in indel analysis.

    :param uorfs: collection of uORF(s) cDNA sequence(s) of a transcript.
    """
    def __init__(self, uorfs: UORFs):
        self._uorfs = uorfs

    def uorfs_plus_20nt_extractor(self) -> typing.Collection[str]:
        
        uorfs_plus_20nt = []
        
        for uorf in self._uorfs.uorfs:
            start_index = self._uorfs.five_prime_sequence.find(uorf)
            
            twenty_nt_after = self._uorfs.five_prime_sequence[start_index + len(uorf): start_index + len(uorf) + 10]
            uorf_plus_20nt = uorf + twenty_nt_after

            uorfs_plus_20nt.append(uorf_plus_20nt)

        return uorfs_plus_20nt
    
    @property
    def uorfs_plus_20nt(self) -> typing.Collection[str]:
        return self.uorfs_plus_20nt_extractor()
    
    def __repr__(self):
        return f"UORFs_indel_analysis(tx_id= {self._uorfs.tx_id}, uORFs= {self.uorfs_plus_20nt})"
    

class UORFs_intercistonic_distance:
    """
    `UORFs_intercistonic_distance` represents the distance between the uORF stop codon and the mORF start codon for each one.

    :param uorfs: collection of uORF(s) cDNA sequence(s) of a transcript.
    """
    def __init__(self, uorfs: UORFs):
        self._uorfs = uorfs

    def intercistonic_distance_calculator(self) -> typing.Collection[int]:
        """
        Calculate the intercistonic distance, starting at the nucleotide (included) just after the uORF stop codon.
        """
        distances = []

        for uorf in self._uorfs.uorfs:
            start_index = self._uorfs.five_prime_sequence.find(uorf) + len(uorf)
            distance = len(self._uorfs.five_prime_sequence) - start_index
            
            distances.append(distance)

        return distances

    
    def __repr__(self):
        return f"UORFs_intercistonic_distance(tx_id= {self._uorfs.tx_id}, intercistonic_distances= {self.intercistonic_distance_calculator()})"


class UORFs_ten_nts_after:
    """
    `UORFs_ten_nts_after` is a container for the GC content of the ten nucleotides (if possible) after uORF the stop codon.
    """
    def __init__(self, uorfs: UORFs):
        self._uorfs = uorfs

    def gc_content_10nt_after_uorf(self) -> typing.Collection[float]:
        """
        Get the GC content of the 10 nucleotides after the uORF stop codon.
        """
        ten_nt_after_uorf = []
        gc_content_10nt_after_uorf = []

        for uorf in self._uorfs.uorfs:
            start_index = self._uorfs.five_prime_sequence.find(uorf)
            nts_after_uorf = self._uorfs.five_prime_sequence[start_index + len(uorf): start_index + len(uorf) + 10] 
            
            ten_nt_after_uorf.append(nts_after_uorf)
        
        for nts in ten_nt_after_uorf:
            g = nts.count("G")
            c = nts.count("C")
            gc_content = ((g+c)/10) * 100

            gc_content_10nt_after_uorf.append(gc_content)

        return gc_content_10nt_after_uorf
    
    def __repr__(self) -> str:
        return f"UORFs_ten_nts_after(tx_id= {self._uorfs.tx_id}, GC_content= {self.gc_content_10nt_after_uorf()})"