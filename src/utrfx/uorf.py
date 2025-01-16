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
        bases = len(five_sequence) - uorf.uorf.end

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
    if uorf.uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")
    
    if bases > (len(five_sequence) - uorf.uorf.end):
        bases = len(five_sequence) - uorf.uorf.end

    return five_sequence[uorf.uorf.start: uorf.uorf.end + bases]


def intercistronic_distance(five_sequence: str, uorf: UORFCoordinates) -> int:
    """
    Calculate the intercistronic distance, defined as the distance from the uORF stop codon to the mORF start codon.

    See here: Silva, J., Fernandes, R., Romão, L. (2019). Translational Regulation by Upstream Open Reading Frames and Human Diseases. 
    In: Romão, L. (eds) The mRNA Metabolism in Human Disease. Advances in Experimental Medicine and Biology, vol 1157. 
    Springer, Cham. https://doi.org/10.1007/978-3-030-19966-1_5
    """
    if uorf.uorf.end > len(five_sequence):
        raise ValueError("uORF overlaps with the mORF")

    return len(five_sequence) - uorf.uorf.end


def cap_five_to_uorf_distance(uorf: UORFCoordinates) -> int:
    """
    Calculate the number of bases between the 5' cap and the uORF start codon.
    """
    return uorf.uorf.start


def kozak_sequence_strength(five_sequence: str, uorf: UORFCoordinates) -> int:
    """
    Indicate the difference between a given Kozak sequence and the consensus sequence, based on two residues.

    The Kozak consensus sequence is defined as CCRCCAUGG, with a purine base in position -3 and a guanine base in position +4 
    as most important for initiation. Initiation sequence contexts are frequently classified as strong (both critical residues match the consensus sequence),
    as adequate/intermediate (either residue -3 or +4 matches) or as weak (neither residue matches).

    :returns: integer being: 2 (strong), 1 (adequate) or 0 (weak).

    See here: Wethmar K, Smink JJ, Leutz A. Upstream open reading frames: molecular switches in (patho)physiology. 
    Bioessays. 2010 Oct;32(10):885-93. doi: 10.1002/bies.201000037. Epub 2010 Aug 19. PMID: 20726009; PMCID: PMC3045505.
    """
    purines = ["A", "G"]
    minus_three_residue = five_sequence[uorf.uorf.start - 3]
    plus_four_residue = five_sequence[uorf.uorf.start + 3]

    if minus_three_residue in purines and plus_four_residue == "G":
        return 2
    elif minus_three_residue in purines and plus_four_residue != "G":
        return 1
    elif minus_three_residue not in purines and plus_four_residue == "G":
        return 1
    else:
        return 0