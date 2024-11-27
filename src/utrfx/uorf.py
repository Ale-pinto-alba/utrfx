import typing

from Bio import SeqIO

class UORFsProcessor:
    """
    `UORFsProcessor` takes a FASTA file with the cDNA sequences of a transcript and its parts to extract different
    features of the uORFs available.

    :param fpath: FASTA file path.
    :ivar seq_records: list of every record of the class 'Bio.SeqRecord.SeqRecord' in the file.
    :ivar tx_record: record (class 'Bio.SeqRecord.SeqRecord') of the transcript.
    :ivar tx_seq: string of the complete transcript nucleotide sequence.
    :ivar tx_id: GENCODE transcript identifier.
    :ivar five_utr_seq: string of the 5'UTR nucleotide sequence.
    :ivar uorfs: list of uORFs (if any).
    :ivar uorfs_with_20nt_more: list of uORFs (if any) with the corresponding 20 nucleotides downstream (if possible).
    """
    def __init__(self, fpath: str):
        self._fpath = fpath
        self._seq_records = self._parse_fasta()
        self._tx_record = self._get_tx_record()
        self._tx_seq = str(self._tx_record.seq).strip()
        self._tx_id = self._tx_record.id
        self._five_utr_seq = self._get_five_utr_sequence()
        self._uorfs = self._uorf_extractor()
        self._uorfs_with_20nt_more = self._uorfs_plus_20nt_extractor()

    @property
    def tx_id(self) -> str:
        return self._tx_id
    
    @property 
    def tx_sequence(self) -> str:
        return self._tx_seq
    
    @property
    def five_utr_sequence(self) -> str:
        return self._five_utr_seq
    
    @property
    def uorfs(self) -> typing.List[str]:
        return self._uorfs
    
    @property
    def uorfs_with_20nt_more(self) -> typing.List[str]:
        return self._uorfs_with_20nt_more

    def five_utr_lenght(self) -> int:
        return len(self._five_utr_seq)
    
    def _parse_fasta(self) -> typing.List:
        seq_records = list(SeqIO.parse(self._fpath, "fasta"))
        if not seq_records:
            raise ValueError("Empty FASTA file or no records.")
        return seq_records

    def _get_tx_record(self):
        for seq_record in self._seq_records:
            if "cdna" in seq_record.description:
                tx_cdna = seq_record    
                if not tx_cdna:
                    raise ValueError("No transcript cDNA in the FASTA file.")
                return tx_cdna

    def _get_five_utr_sequence(self) -> str:
        for seq_record in self._seq_records:
            if "utr5" in seq_record.description:
                five_utr_seq = str(seq_record.seq).strip()   
                assert isinstance(five_utr_seq, str)
                if not five_utr_seq:
                    raise ValueError("No 5'UTR region in the FASTA file.")
                return five_utr_seq
            
    def _uorf_extractor(self) -> typing.List[str]:
        """
        Take the nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.
        """
        five_sequence = self._five_utr_seq
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
            # It resumes the search a codon after the start codon of the current uORF, but should I resume it after the stop
            # codon or I don´t look for overlapping uORFs?
            
        return uorfs
    
    def _uorfs_plus_20nt_extractor(self) -> typing.List[str]:
        """
        Search the uORFs in their 5'UTR region to extract the uORFs with the 20 nucleotides (if possible) after the
        corresponding stop codon for indel analysis.
        """
        uorfs_plus_20nt = []

        for uorf in self._uorfs:
            start_index = self._five_utr_seq.find(uorf)
            twenty_nts_after_uorf = self._five_utr_seq[start_index + len(uorf): start_index + len(uorf) + 20] 
            
            uorf_plus_20nt = uorf + twenty_nts_after_uorf
            uorfs_plus_20nt.append(uorf_plus_20nt)
            assert uorf in uorf_plus_20nt, "uORF lost."

        return uorfs_plus_20nt

    def number_of_uorfs(self) -> int:
        return len(self._uorfs)
    
    def uorfs_lengths(self) -> typing.List[int]:
        return [len(uorf) for uorf in self._uorfs]
    
    def gc_content(self) -> typing.List[float]:
        """
        Get the GC content of each uORF.
        """
        gc_content_per_uorf = []

        for uorf in self._uorfs:
            total = len(uorf)
            g = uorf.count("G")
            c = uorf.count("C")
            
            gc_content = ((g+c) / total) * 100
            gc_content_per_uorf.append(gc_content)

        return gc_content_per_uorf
    
    def intercistonic_distance(self) -> typing.List[int]:
        """
        Calculate the intercistonic distance (distance between the uORF stop codon and the mORF start codon).
        """
        distances = []

        for uorf in self._uorfs:
            start_index = self._five_utr_seq.find(uorf) + len(uorf)
            distance = len(self._five_utr_seq) - start_index
            
            distances.append(distance)
        return distances
    
    def gc_content_10nt_after_uorf(self) -> typing.List[float]:
        """
        Get the GC content of the 10 nucleotides (if possible) after the uORF stop codon.
        """
        ten_nt_after_uorf = []
        gc_content_10nt_after_uorf = []

        for uorf in self._uorfs:
            start_index = self._five_utr_seq.find(uorf)
            nts_after_uorf = self._five_utr_seq[start_index + len(uorf): start_index + len(uorf) + 10] 
            
            ten_nt_after_uorf.append(nts_after_uorf)
        
        for nts in ten_nt_after_uorf:
            g = nts.count("G")
            c = nts.count("C")
            gc_content = ((g+c)/10) * 100

            gc_content_10nt_after_uorf.append(gc_content)
        return gc_content_10nt_after_uorf

    def __repr__(self) -> str:
        return f"UORFsProcessor(tx_id= {self.tx_id}, uORFs= {self.uorfs})"