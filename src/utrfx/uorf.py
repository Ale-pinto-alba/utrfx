import abc
import typing

import requests

from utrfx.genome import Region, GenomeBuild, GenomicRegion, Strand, GRCh38
from utrfx.model import FiveUTR, FivePrimeSequence, UORF, Transcript

def download_fasta_from_ensembl(transcript_id: str, output_file: str) -> typing.IO:
    """
    Download a FASTA file containing the cDNA sequence for a given transcript from Ensembl's REST API.

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`
    :param output_file: the path to save the downloaded FASTA file.
    """
    base_url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?content-type=text/x-fasta"

    response = requests.get(base_url)

    if response.status_code == 200:
        with open(output_file, 'w') as f:
            f.write(response.text)
        print(f"FASTA file successfully downloaded and saved to {output_file}")
    else:
        raise Exception(f"Error: Unable to download FASTA. HTTP Status Code: {response.status_code}")
    

def get_five_prime_sequence(transcript_sequence: str, five_utrs: FiveUTR) -> str:
    """
    Return the 5'UTR cDNA sequence of a given transcript nucleotide sequence.

    :param transcript_sequence: FASTA file with the transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    """
    with open(transcript_sequence) as f:
        lines = f.readlines()
        sequence = "".join(line.strip() for line in lines[1:])
        return sequence[:five_utrs.__len__()]


def uorf_extractor(five_utr: FivePrimeSequence) -> typing.Collection[UORF]:
    """
    Take the nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.
    """
    five_sequence = five_utr._five_prime_sequence
    uorfs = []
    start_position = 0

    while start_position < len(five_sequence) - 3:
        start_index = five_sequence.find("ATG", start_position)
        if start_index < 0:
            break  # No more uORF in the remaining sequence

        found_stop = False
        for i in range(start_index, len(five_sequence) - 2, 3):  
            codon = five_sequence[i:i + 3]
            if codon in ["TAA", "TAG", "TGA"]:
                stop_index = i + 3
                found_stop = True
                break

        if found_stop:
            uorfs.append(UORF(
                five_prime=five_utr,
                uorf=Region(start=start_index, end=stop_index)
            ))
            start_position = stop_index  
        else:
            start_position = start_index + 3  

    return uorfs


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
    
class GcContentCalculator:

    def __init__(
        self,
        uorf: UORF,
    ):
        self._uorf = uorf

    def gc_content_n_bases_downstream(self, bases: int) -> float:
        """
        """
        if bases > (len(self._uorf._five_prime.five_prime_sequence) - self._uorf._uorf.end):
            bases = len(self._uorf._five_prime.five_prime_sequence) - self._uorf._uorf.end
            print(f"Limit exceed, calculating the GC content of the {total} nucleotides after the uORF")

        region_downstream = Region(start= self._uorf._uorf.end, end= self._uorf._uorf.end + bases)
        total = region_downstream.__len__()
        if total == 0:
            return "GC content equals 0"
        else:
            g = self._uorf._five_prime.five_prime_sequence[region_downstream.start: region_downstream.end].count("G")
            c = self._uorf._five_prime.five_prime_sequence[region_downstream.start: region_downstream.end].count("C")
            gc_content = ((g+c)/total) * 100

            return gc_content


def uorfs_plus_n_nts_downstream_extractor(uorf: UORF, bases: int) -> str:
    """ 
    Get the uORF plus the `n` nucleotides after the uORF stop codon (if possible) for indel analysis.
    """
    if bases > (len(uorf._five_prime.five_prime_sequence) - uorf._uorf.end):
            bases = len(uorf._five_prime.five_prime_sequence) - uorf._uorf.end
            print(f"Limit exceed, taking the uORF plus the {bases} nucleotides downstream")

    region_downstream = Region(start= uorf._uorf.end, end= uorf._uorf.end + bases)

    return uorf.uorf_sequence + uorf._five_prime.five_prime_sequence[region_downstream.start: region_downstream.end]


class KozakSequenceCalculator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def compute(self) -> float:
        pass


class FooKozakCalculator(KozakSequenceCalculator):
    
    def compute(self) -> float:
        return 5.