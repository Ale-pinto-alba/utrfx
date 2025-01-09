import typing

import requests

from utrfx.genome import Region
from utrfx.model import FiveUTRCoordinates, UORFCoordinates


def download_fasta_from_ensembl(transcript_id: str, timeout: float = 30.,) -> str:
    """
    Download a FASTA file containing the cDNA sequence (spliced mRNA) for a given transcript from Ensembl's REST API and
    return only the nucleotide sequence (without the FASTA header).

    If the transcript is in the negative strand, it will be given the reverse complement of the spliced mRNA found in Ensembl's REST API,
    which is it located in the positive strand.         

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`
    """
    base_url_cdna = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna&content-type=text/x-fasta"
    
    response = requests.get(base_url_cdna, timeout=timeout)

    if response.status_code == 200:
        lines = response.text.splitlines()

        cdna_sequence = ''.join(lines[1:])

        base_url_strand = f"https://rest.ensembl.org/lookup/id/{transcript_id}?content-type=application/json"

        response_strand = requests.get(base_url_strand, timeout=timeout)
        data = response_strand.json()

        strand = data.get("strand", None)
        
        if strand > 0: 
            return cdna_sequence
        elif strand < 0:
            return cdna_sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1]
        else:
            raise Exception("Not strand found.")
    else:
        response.raise_for_status()


def get_five_prime_sequence(cdna_sequence: str, five_utrs: FiveUTRCoordinates) -> str:
    """
    Return the 5'UTR cDNA sequence of a given transcript nucleotide sequence.

    :param transcript_sequence: transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    """
    return cdna_sequence[:len(five_utrs)]


def uorf_extractor(five_utr: FiveUTRCoordinates, five_sequence: str) -> typing.Collection[UORFCoordinates]:
    """
    Take the cDNA nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.

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