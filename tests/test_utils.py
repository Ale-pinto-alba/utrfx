import pytest

from utrfx.util import download_fasta_from_ensembl

@pytest.mark.online
def test_download_fasta_from_ensembl():
    fasta = download_fasta_from_ensembl("ENST00000381418")

    assert fasta.startswith("AGTT")

    assert len(fasta) > 0