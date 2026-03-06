import typing

import requests
import pandas as pd

from utrfx.genome import *
from utrfx.model import *
from utrfx.gtf_io import GTFio
from utrfx.uorf import *
from utrfx.variant_util import *
from utrfx.rna_util import *


def fetch_cdna_from_ensembl(transcript_id: str, timeout: float = 30.,) -> str:
    """
    Download cDNA sequence (spliced mRNA) for a given transcript from Ensembl's REST API.

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`.
    :param timeout: timeout for the API request in seconds (default: 30 seconds).
    :returns: a `str` with the cDNA sequence.
    """
    base_url_cdna = f"https://rest.ensembl.org/sequence/id/{transcript_id}?type=cdna&content-type=text/x-fasta"
    
    response = requests.get(base_url_cdna, timeout=timeout)

    if response.status_code == 200:
        lines = response.text.splitlines()
        return ''.join(lines[1:])
    else:
        response.raise_for_status()


def get_five_prime_sequence(cdna_sequence: str, five_utrs: FiveUTRCoordinates) -> str:
    """
    Return the 5'UTR cDNA sequence of a given transcript nucleotide sequence (spliced mRNA).

    :param transcript_sequence: transcript nucleotide sequence.
    :param five_utrs: 5'UTR Genomic Region(s).
    :returns: a `str` with the 5'UTR cDNA sequence.
    """
    return cdna_sequence[:len(five_utrs)]


def uorf_extractor(five_utr: FiveUTRCoordinates, five_sequence: str) -> typing.Collection[UORFCoordinates]:
    """
    Take the cDNA nucleotide sequence of a transcript 5'UTR region to extract the uORFs sequences available.

    :param five_utr: list of Genomic Regions corresponding to the 5'UTRs regions.
    :param five_sequence: 5'UTR cDNA sequence.
    :returns: a list of UORFCoordinates objects, each containing the uORF region and whether it is an overlapping uORF (ouORF).
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
                uorf=Region(start=start_index, end=stop_index),
                ouorf= False,
            ))
            start_position = stop_index  
        else:
            uorfs.append(UORFCoordinates(
                five_utr=five_utr,
                uorf=Region(start=start_index, end=len(five_sequence)),
                ouorf= True,
            ))
            start_position = start_index + 1  

    return uorfs

def sum_all_features(row: pd.Series) -> float:
    """
    Sum all numeric values in a row (ignoring non-numeric columns).
    
    :param row: A row from a pandas DataFrame.
    :returns: The sum of all numeric values in the row.
    """
    numeric_row = pd.to_numeric(row, errors='coerce')
    return numeric_row.sum()

def obtain_all_features_as_list(
    wt_five_utr_cdna_sequence: str,
    variant_five_utr_cdna_sequence: str,
) -> typing.List:
    """
    :param wt_five_utr_cdna_sequence: The wild-type 5' UTR cDNA sequence.
    :param variant_five_utr_cdna_sequence: The variant 5' UTR cDNA sequence.
    :returns: A list of all features.
    """

    wt_fc = RNA.fold_compound(wt_five_utr_cdna_sequence)
    wt_structure, wt_mfe = wt_fc.mfe()
    wt_fc.pf()
    wt_canonical_diversity = wt_fc.mean_bp_distance()
    wt_mfe_freq = wt_fc.pr_structure(wt_structure)
    wt_lbox, wt_ubox = generate_probs(wt_five_utr_cdna_sequence)

    variant_fc = RNA.fold_compound(variant_five_utr_cdna_sequence)
    variant_structure, variant_mfe = variant_fc.mfe()
    variant_fc.pf()
    variant_canonical_diversity = variant_fc.mean_bp_distance()
    variant_mfe_freq = variant_fc.pr_structure(variant_structure)
    variant_lbox, variant_ubox = generate_probs(variant_five_utr_cdna_sequence)