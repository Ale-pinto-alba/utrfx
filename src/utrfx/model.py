import typing
    
from .genome import GenomicRegion, Region


class FiveUTR:
    """
    `FiveUTR` is a container for 5'UTR Genomic Regions.

    :param regions: list of Genomic Regions corresponding to the 5'UTRs regions.
    """
    def __init__(
        self,
        regions: typing.Collection[GenomicRegion],
    ):
        self._regions = regions

    @property
    def regions(self) -> typing.Collection[GenomicRegion]:
        return self._regions

    def __len__(self) -> int:
        return 0 if len(self._regions) == 0 else sum(len(region) for region in self._regions)
    
    def __repr__(self):
        regions_info = ", ".join([f"({region._contig.ucsc_name}, {region.start}, {region.end}, {region.strand})" for region in self._regions])
        return f"FiveUTR(regions={len(self._regions)} regions: {regions_info})"
    

class Transcript:
    """
    `Transcript` represents the 5'UTR Genomic Region(s) of a transcript.

    :param tx_id: GENCODE transcript identifier, e.g. `ENST00000381418`
    :param five_utr: 5'UTR Genomic Region(s). 
    """
    def __init__(
        self, 
        tx_id: str, 
        five_utr: FiveUTR
    ):
        self._tx_id = tx_id
        self._five_utr = five_utr

    @property
    def tx_id(self) -> str:
        return self._tx_id
    
    @property
    def five_utr(self) -> FiveUTR:
        return self._five_utr

    def __repr__(self):
        return f"Transcript(tx_id={self._tx_id}, five_utr={repr(self._five_utr)})"


class FivePrimeSequence:
    """
    `FivePrimeSequence` represents the 5'UTR cDNA sequence of a transcript.

    :param transcript: transcript with its GENCODE ID and 5'UTR Genomic Region(s).
    :param five_prime_sequence: 5'UTR nucleotide sequence.
    """
    def __init__(
        self,
        transcript: Transcript,
        five_prime_sequence: str,
    ):
        self._transcript = transcript
        self._five_prime_sequence = five_prime_sequence.upper()

    @property
    def transcript(self) -> typing.Sequence[GenomicRegion]:
        return self._transcript

    @property
    def five_prime_sequence(self) -> str:
            return self._five_prime_sequence
        
    # @staticmethod
    # def from_fasta(tx: Transcript, fpath: str) -> "FivePrimeSequence":
    #     valid_fasta_extensions = (".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn")

    #     if fpath.endswith(valid_fasta_extensions):
    #         with open(fpath) as f:
    #             lines = f.readlines()
    #             return FivePrimeSequence(
    #                 tx= Transcript(tx_id= tx.tx_id, five_utr= tx.five_utr),
    #                 five_prime_sequence= "".join(line.strip() for line in lines[1:]),
    #             )
    #     else:
    #         raise ValueError("Unexpected extension.")

    def __len__(self) -> int:
        return len(self._five_prime_sequence)
    
    def __eq__(self, other):
        return (isinstance(other, FivePrimeSequence)
                and self._transcript == other._transcript
                and self._five_prime_sequence == other._five_prime_sequence)

    def __repr__(self) -> str:
        return f"FivePrimeSequence(tx_id= {self._transcript.tx_id}, 5'UTRs= {self._transcript.five_utr}, 5'UTR_sequence= {self.five_prime_sequence})"
    

class UORF:
    """
    `UORF` represents an uORF cDNA sequence of a transcript.

    The UORF is upstream of the mORF and they do *not* overlap.

    :param five_prime: transcript with its GENCODE ID, 5'UTR Genomic Region(s) and cDNA sequence.
    :param uorf: uORF region marked by its start and end nucleotide.
    """
    def __init__(
        self,
        five_prime: FivePrimeSequence,
        uorf: Region,
    ):
        self._five_prime = five_prime
        self._uorf = uorf

    @property
    def uorf_sequence(self) -> str:
        return self._five_prime.five_prime_sequence[self._uorf.start: self._uorf.end]

    def intercistonic_distance(self) -> int:
        """
        Calculate the intercistonic distance, which is the number of bases located between the uORF stop codon 
        and the mORF start codon, starting at the nucleotide (included) just after the uORF stop codon.
        """
        intercistonic_distance = len(self._five_prime.five_prime_sequence) - self._uorf.end

        return intercistonic_distance

    def __len__(self) -> int:
        return len(self.uorf_sequence)  

    def __eq__(self, other):
        return (isinstance(other, UORF)
                and self._five_prime== other._five_prime
                and self._uorf == other._uorf)
    
    def __repr__(self) -> str:
        return f"UORF(tx_id= {self._five_prime._transcript.tx_id}, 5'UTRs= {self._five_prime._transcript.five_utr}, 5'UTR_sequence= {self._five_prime.five_prime_sequence}, uORF= {self._uorf})"
