import abc
import typing

from utrfx.genome import Region
from utrfx.model import FiveUTRCoordinates, UORFCoordinates


def get_five_prime_sequence(transcript_sequence: str, five_utrs: FiveUTRCoordinates) -> str:
    """
    Return the 5'UTR cDNA sequence of a given transcript nucleotide sequence.

    :param transcript_sequence: transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    """
    return transcript_sequence[:len(five_utrs)]


def uorf_extractor(five_utr: FiveUTRCoordinates, five_sequence: str) -> typing.Collection[UORFCoordinates]:
    """
    Take the nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.

    :param five_utr: list of Genomic Regions corresponding to the 5'UTRs regions.
    :param five_sequence: 5'UTR cDNA sequence.
    """
    uorfs = []
    start_position = 0

    while start_position < len(five_sequence) - 2:
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
            uorfs.append(UORFCoordinates(
                five_utr=five_utr,
                uorf=Region(start=start_index, end=stop_index)
            ))
            start_position = stop_index  
        else:
            start_position = start_index + 3  

    return uorfs


def gc_content(five_sequence: str, uorf: UORFCoordinates) -> float:
    """
    Get the GC content of an uORF.

    If the number of bases chosen ends beyond the 5'UTR end limit (overlapping with the mORF), the number of bases taken will be clipped to
    the length between the uORF stop codon and the mORF start codon.

    :returns: the GC content as a float in range [0, 1].
    """
    total = uorf.uorf.end - uorf.uorf.start

    if uorf.uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")

    if total == 0:
        return 0
    else:
        region = five_sequence[uorf.uorf.start: uorf.uorf.end]
        g = region.count("G")
        c = region.count("C")
                
        return (g+c) / total
    

def gc_content_n_bases_downstream(five_sequence: str, uorf: UORFCoordinates, bases: int) -> float:
    """
    Get the GC content of a `n` bases downstream an uORF.

    If the number of bases chosen ends beyond the 5'UTR end limit (overlapping with the mORF), the number of bases taken will be clipped to
    the length between the uORF stop codon and the mORF start codon.

    :param bases: a non-negative number of downstream bases

    :returns: the GC content as a float in range [0, 1].
    """
    if uorf.uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")
    
    assert bases >= 0
    
    if bases > (len(five_sequence) - uorf.uorf.end):
        bases = (len(five_sequence) - uorf.uorf.end)

    total = uorf.uorf.end + bases - uorf.uorf.end
    if total == 0:
        return 0
    else:
        region = five_sequence[uorf.uorf.end: uorf.uorf.end + bases]
        g = region.count("G")
        c = region.count("C")
        
        return (g+c)/total



def uorfs_plus_n_nts_downstream_extractor(five_sequence: str, uorf: UORFCoordinates, bases: int) -> str:
    """ 
    Get the uORF plus the `n` nucleotides downstream of the uORF stop codon (if possible) for indel analysis.

    If the number of bases chosen ends beyond the 5'UTR end limit (overlapping with the mORF), the number of bases taken will be clipped to
    the length between the uORF stop codon and the mORF start codon.
    """
    if uorf._uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")
    
    if bases > (len(five_sequence) - uorf._uorf.end):
        bases = len(five_sequence) - uorf._uorf.end

    return five_sequence[uorf.uorf.start: uorf.uorf.end + bases]


def intercistonic_distance(five_sequence: str, uorf: UORFCoordinates) -> int:
    """
    Calculate the intercistonic distance, which is the number of bases located between the uORF stop codon 
    and the mORF start codon.
    """
    if uorf._uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")

    return len(five_sequence) - uorf.uorf.end


class KozakSequenceCalculator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def compute(self) -> float:
        pass


class FooKozakCalculator(KozakSequenceCalculator):
    
    def compute(self) -> float:
        return 5.
