import abc
import typing

import requests

from utrfx.genome import Region
from utrfx.model import FiveUTR, FivePrimeSequence, UORF

def download_fasta_from_ensembl(transcript_id: str) -> str:
    """
    Download a FASTA file containing the cDNA sequence for a given transcript from Ensembl's REST API.

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`
    """
    base_url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?content-type=text/x-fasta"
    try:
        response = requests.get(base_url, timeout= 30)
    except requests.exceptions.Timeout:
        print("Timed out")

    if response.status_code == 200:
        return response.text[1:]
    else:
        raise Exception(f"Error: Unable to download FASTA. HTTP Status Code: {response.status_code}")
    

def get_five_prime_sequence(transcript_sequence: str, five_utrs: FiveUTR) -> str:
    """
    Return the 5'UTR cDNA sequence of a given transcript nucleotide sequence.

    :param transcript_sequence: transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    """
    return transcript_sequence[:len(five_utrs)]


def uorf_extractor(five_utrs: FiveUTR, five_utr_seq: str) -> typing.Collection[UORF]:
    """
    Take the nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.

    :param five_utr: list of Genomic Regions corresponding to the 5'UTRs regions.
    :param seq: 5'UTR cDNA sequence.
    """
    five_sequence = five_utr_seq
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
                five_utr=five_utrs,
                uorf=Region(start=start_index, end=stop_index)
            ))
            start_position = stop_index  
        else:
            start_position = start_index + 3  

    return uorfs


def gc_content(five_utr_seq: str, uorf: UORF) -> float:
    """
    Get the GC content of an uORF.
    """
    total = uorf.__len__()

    if total == 0:
        raise ValueError("GC content is 0")
    else:
        g = five_utr_seq[uorf._uorf.start: uorf._uorf.end].count("G")
        c = five_utr_seq[uorf._uorf.start: uorf._uorf.end].count("C")
                
        gc_content = ((g+c) / total) * 100

        return gc_content
    

def gc_content_n_bases_downstream(five_utr_seq: str, uorf: UORF, bases: int) -> float:
    """
    Get the GC content of a `n` bases downstream an uORF.
    """
    if bases > (len(five_utr_seq) - uorf._uorf.end):
        bases = (len(five_utr_seq) - uorf._uorf.end)
        print(f"Limit exceed, calculating the GC content of the {total} nucleotides after the uORF")

    region_downstream = Region(start= uorf._uorf.end, end= uorf._uorf.end + bases)
    total = region_downstream.__len__()
    if total == 0:
        return "GC content equals 0"
    else:
        g = five_utr_seq[region_downstream.start: region_downstream.end].count("G")
        c = five_utr_seq[region_downstream.start: region_downstream.end].count("C")
        gc_content = ((g+c)/total) * 100

        return gc_content


def uorfs_plus_n_nts_downstream_extractor(uorf: UORF, five_utr_seq: str, bases: int) -> str:
    """ 
    Get the uORF plus the `n` nucleotides after the uORF stop codon (if possible) for indel analysis.
    """
    if bases > (len(five_utr_seq) - uorf._uorf.end):
            bases = len(five_utr_seq) - uorf._uorf.end
            print(f"Limit exceed, taking the uORF plus the {bases} nucleotides downstream")

    region_downstream = Region(start= uorf._uorf.end, end= uorf._uorf.end + bases)

    return five_utr_seq[uorf._uorf.start:uorf._uorf.end] + five_utr_seq[region_downstream.start: region_downstream.end]


def intercistonic_distance(five_utr_seq: str, uorf: UORF) -> int:
    """
    Calculate the intercistonic distance, which is the number of bases located between the uORF stop codon 
    and the mORF start codon, starting at the nucleotide (included) just after the uORF stop codon.
    """
    intercistonic_distance = len(five_utr_seq) - uorf._uorf.end

    return intercistonic_distance


class KozakSequenceCalculator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def compute(self) -> float:
        pass


class FooKozakCalculator(KozakSequenceCalculator):
    
    def compute(self) -> float:
        return 5.