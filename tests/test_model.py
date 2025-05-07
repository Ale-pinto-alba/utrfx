import os
import pytest

from utrfx.model import TxperGene

@pytest.fixture(scope="module")
def json_fpath(fpath_data_dir: str) -> str:
    return os.path.join(fpath_data_dir, "Ensembl_transcript_per_gene_dictionary.json")

@pytest.fixture(scope="module")
def tx_per_gene_class(json_fpath: str) -> TxperGene:
    return TxperGene(fpath=json_fpath)


@pytest.mark.parametrize(
    "gene_symbol, expected",
    [
        ("ADSL", "ENST00000623063.3"),
        ("HR", "ENST00000381418.9"),
        ("LPL", "ENST00000650287.1"),
        ("NOTEXIST", None)
    ]
)
def test_gene_symbol(
    gene_symbol: str,
    expected: str,
    tx_per_gene_class: TxperGene,
):
    actual = tx_per_gene_class.get_transcript_id(gene_symbol=gene_symbol)
    assert actual == expected