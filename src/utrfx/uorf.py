import abc

class FivePrimeSequence:
    """
    `FivePrimeSequence` represents the 5'UTR cDNA sequence of a transcript.

    :param tx_id: GENCODE transcript identifier e.g. `ENST00000381418_9`
    :param five_prime_sequence: 5'UTR nucleotide sequence.
    """
    def __init__(self, tx_id: str, five_prime_sequence: str):
        self._tx_id = tx_id
        self._five_prime_sequence = five_prime_sequence

    @property
    def tx_id(self) -> str:
        return self._tx_id
    
    @property
    def five_prime_sequence(self) -> str:
            return self._five_prime_sequence
        
    @staticmethod
    def from_fasta(tx_id: str, fpath: str) -> "FivePrimeSequence":
        valid_fasta_extensions = (".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn")

        if fpath.endswith(valid_fasta_extensions):
            with open(fpath) as f:
                lines = f.readlines()
                return FivePrimeSequence(tx_id=tx_id, five_prime_sequence= "".join(line.strip() for line in lines[1:]))
        else:
            raise ValueError("Unexpected extension.")

    def __len__(self) -> int:
        return len(self._five_prime_sequence)
    
    def __eq__(self, other):
        return (isinstance(other, FivePrimeSequence)
                and self.tx_id == other.tx_id
                and self.five_prime_sequence == other.five_prime_sequence)

    def __repr__(self) -> str:
        return f"FivePrimeSequence(tx_id= {self.tx_id}, 5'UTR= {self.five_prime_sequence})"


class UORF:
    """
    `UORF` represents an uORF cDNA sequence of a transcript. 

    :param tx_id: GENCODE transcript identifier e.g. `ENST00000381418_9`
    :param uorf: uORF nucleotide sequence.
    """
    def __init__(self, tx_id: str, uorf: str):
        self._tx_id = tx_id
        self._uorf = uorf

    @property
    def tx_id(self) -> str:
        return self._tx_id
    
    @property
    def uorf_sequence(self) -> str:
        return self._uorf   
    
    def is_in_five_utr(self, five_utr: FivePrimeSequence) -> bool:
        """
        Check if the uORF belongs to a given 5'UTR sequence.
        """
        if self._uorf in five_utr.five_prime_sequence:
            return True
        else:
            return False
    
    def __len__(self) -> int:
        return len(self._uorf)  

    def __eq__(self, other):
        return (isinstance(other, UORF)
                and self.tx_id == other.tx_id
                and self.uorf_sequence == other.uorf_sequence)
    
    def __repr__(self) -> str:
        return f"UORF(tx_id= {self.tx_id}, uORF= {self.uorf_sequence})"


def extract_next_uorf(five_prime_sequence: FivePrimeSequence) -> UORF:
    """
    Extract the next uORF from the 5'UTR sequence and update it.
    Returns the uORF as a UORF instance or None if no more are found.
    """
    five_sequence = five_prime_sequence.five_prime_sequence.upper()
    stop_codons = {"TAA", "TAG", "TGA"} 

    if "ATG" not in five_sequence:
        return None

    start_index = five_sequence.find("ATG")
    codon_list = []

    for i in range(start_index, len(five_sequence) - len(five_sequence) % 3, 3):
        codon = five_sequence[i:i + 3]

        if len(codon) < 3:
            break

        codon_list.append(codon)

        if codon in stop_codons:
            five_sequence = five_sequence[i + 3:]
            return UORF(tx_id= five_prime_sequence.tx_id, uorf= "".join(codon_list))
    

def gc_content(uorf: UORF) -> float:
    """
    Get the GC content of an uORF.
    """
    total = len(uorf.uorf_sequence)

    if total == 0:
        raise ValueError("Cannot calculate GC content from a 0 nt uORF.")
    else:
        g = uorf.uorf_sequence.count("G")
        c = uorf.uorf_sequence.count("C")
                
        gc_content = ((g+c) / total) * 100

        return gc_content
    

def intercistonic_distance(five_utr: FivePrimeSequence, uorf: UORF) -> float:
    """
    Calculate the intercistonic distance, which is the distance between the uORF stop codon and the mORF start codon, 
    starting at the nucleotide (included) just after the uORF stop codon.
    """
    if uorf.is_in_five_utr(five_utr= five_utr) == False or five_utr.tx_id != uorf.tx_id:
        raise ValueError("The uORF is not in the given 5'UTR region.")
    else:
        start_index = five_utr.five_prime_sequence.find(uorf.uorf_sequence) + len(uorf.uorf_sequence)
        intercistonic_distance = float(len(five_utr.five_prime_sequence) - start_index)

        return intercistonic_distance


def gc_content_of_the_ten_nts_after(five_utr: FivePrimeSequence, uorf: UORF) -> float:
    """
    """
    if uorf.is_in_five_utr(five_utr= five_utr) == False or five_utr.tx_id != uorf.tx_id:
        raise ValueError("The uORF is not in the given 5'UTR region.")
    else:
        start_index = five_utr.five_prime_sequence.find(uorf.uorf_sequence) + len(uorf.uorf_sequence)
        ten_nt_after = five_utr.five_prime_sequence[start_index: start_index + 10]
        
        total = len(ten_nt_after)
        if total == 0:
            raise ValueError("Cannot calculate GC content from a 0 nt region.")
        
        elif total < 10:
            print(f"Calculating the GC content of the {total} nucleotides after the uORF")

            g = ten_nt_after.count("G")
            c = ten_nt_after.count("C")
            gc_content = ((g+c)/total) * 100
            return gc_content
        
        else:

            g = ten_nt_after.count("G")
            c = ten_nt_after.count("C")
            gc_content = ((g+c)/10) * 100
            return gc_content


def uorfs_plus_20nt_extractor(five_utr: FivePrimeSequence, uorf: UORF) -> str:
    """ 
    Get the uORF plus the twenty nucleotides after the uORF stop codon (if possible) for indel analysis.
    """
    if uorf.is_in_five_utr(five_utr= five_utr) == False or five_utr.tx_id != uorf.tx_id:
        raise ValueError("The uORF is not in the given 5'UTR region.")
    else:
        start_index = five_utr.five_prime_sequence.find(uorf.uorf_sequence) + len(uorf.uorf_sequence)
        twenty_nt_after = five_utr.five_prime_sequence[start_index: start_index + 20]

        if twenty_nt_after != 20:
            print(f"Taking the {twenty_nt_after} nucleotides after instead of 20.")
            uorf_plus_nt = uorf.uorf_sequence + twenty_nt_after
            
            return uorf_plus_nt
        else:
            uorf_plus_20nt = uorf.uorf_sequence + twenty_nt_after
        
            return uorf_plus_20nt


class KozakSequenceCalculator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def compute(self) -> float:
        pass


class FooKozakCalculator(KozakSequenceCalculator):
    
    def compute(self) -> float:
        return 5.