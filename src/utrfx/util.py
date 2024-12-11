import requests


def download_fasta_from_ensembl(transcript_id: str, timeout: float = 30.,) -> str:
    """
    Download a FASTA file containing the cDNA sequence for a given transcript from Ensembl's REST API and
    return the only the nucleotide sequence (without the FASTA header).

    :param transcript_id: Ensembl transcript identifier e.g. `ENST00000381418`
    """
    base_url = f"https://rest.ensembl.org/sequence/id/{transcript_id}?content-type=text/x-fasta"
    
    response = requests.get(base_url, timeout=timeout)

    if response.status_code == 200:
        lines = response.text.splitlines()
        return ''.join(lines[1:])
    else:
        response.raise_for_status()