import os

import pytest

from utrfx.uorf import FivePrimeSequence, UORF, extract_next_uorf, uorfs_plus_20nt_extractor, intercistonic_distance, gc_content, gc_content_of_the_ten_nts_after

@pytest.fixture(scope="session")
def fpath_fasta(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data", "Homo_sapiens_ENST00000381418_9_sequence_sample.fa")

@pytest.fixture(scope="session")
def example_five_utr(fpath_fasta: str):
    return FivePrimeSequence.from_fasta(tx_id= "ENST00000381418.9", fpath= fpath_fasta)

def test_five_prime_sequence_length(example_five_utr: FivePrimeSequence):

    assert example_five_utr.__len__() == 623

def test_all_are_nucleotides(example_five_utr: FivePrimeSequence):
    nucleotides = ["A", "T", "C", "G"]

    for nt in example_five_utr.five_prime_sequence.upper():
        
        assert nt in nucleotides

@pytest.fixture
def example_uorf(example_five_utr: FivePrimeSequence):
    return extract_next_uorf(five_prime_sequence= example_five_utr)

def test_uorf(example_uorf: UORF, example_five_utr: FivePrimeSequence):

    assert example_uorf.__len__() == 51

    assert example_uorf.is_in_five_utr(five_utr= example_five_utr) == True

def test_gc_content(example_uorf: UORF):

    assert gc_content(uorf= example_uorf) == 62.745098039215684

def test_intercistonic_distance(example_five_utr: FivePrimeSequence, example_uorf: UORF):

    assert intercistonic_distance(five_utr= example_five_utr, uorf= example_uorf) == 556.0

def test_gc_content_10_nt_after(example_five_utr: FivePrimeSequence, example_uorf: UORF):

    assert gc_content_of_the_ten_nts_after(five_utr= example_five_utr, uorf= example_uorf) == 90.0

def test_uorf_plus_20_nt(example_five_utr: FivePrimeSequence, example_uorf: UORF):

    assert uorfs_plus_20nt_extractor(five_utr= example_five_utr, uorf= example_uorf) == "ATGGCGATCAGAGGTCCTGCTGCGCTCTCCGCCGCGCTCTACCTCCATTAGCCGCGCTGCGCGGTGCTGCG"