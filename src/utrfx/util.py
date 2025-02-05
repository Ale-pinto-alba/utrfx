import typing

import requests

from utrfx.genome import Region
from utrfx.model import FiveUTRCoordinates, UORFCoordinates


def fetch_genomic_sequence_from_ensembl(transcript_id: str, timeout: float = 30.,) -> str:
    """
    Download the genomic sequence for a given transcript from Ensembl's REST API.

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`
    """
    base_url_cdna = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=genomic&content-type=text/x-fasta"
    
    response = requests.get(base_url_cdna, timeout=timeout)

    if response.status_code == 200:
        lines = response.text.splitlines()
        return ''.join(lines[1:])
    else:
        response.raise_for_status()


def get_five_prime_sequence(genomic_sequence: str, five_utrs: FiveUTRCoordinates) -> str:
    """
    Return the 5'UTR DNA sequence of a given transcript nucleotide sequence.

    :param transcript_sequence: transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    """
    five_utrs_tuples = []
    five_utr_sequences = []

    for region in five_utrs.regions:
        five_utrs_tuples.append((region.start, region.end))
    
    five_utrs_tuples.sort()
    previous_end = None
    accumulative_length = 0
    for start,end in five_utrs_tuples:
        five_utr_length = end - start
        if previous_end is None:
            five_utr_sequences.append(genomic_sequence[:five_utr_length])
            previous_end = end
            accumulative_length += five_utr_length

        elif previous_end is not None:
            non_five_utr_region = start - previous_end
            five_utr_sequences.append(genomic_sequence[accumulative_length + non_five_utr_region:accumulative_length + non_five_utr_region + five_utr_length])
            previous_end = end
            accumulative_length += five_utr_length

    return ''.join(five_utr_sequences)


def uorf_extractor(five_utr: FiveUTRCoordinates, five_sequence: str) -> typing.Collection[UORFCoordinates]:
    """
    Take the cDNA nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available
    (not those overlapping with the main ORF).

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