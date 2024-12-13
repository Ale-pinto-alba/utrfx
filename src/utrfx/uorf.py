import abc

from utrfx.model import UORFCoordinates


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
