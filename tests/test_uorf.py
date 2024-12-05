import os

import pytest

from utrfx.uorf import UORFs, UORFs_calculations, UORFs_indel_analysis, UORFs_intercistonic_distance, UORFs_ten_nts_after

@pytest.fixture(scope="session")
def fpath_fasta(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "data", "Homo_sapiens_ENST00000381418_9_sequence_sample.fa")

def test_utr_length(example_uorfs: UORFs):

    assert len(example_uorfs.five_prime_sequence) == 623

def test_uorfs_calculations(example_uorfs: UORFs):

    calculations = UORFs_calculations(uorfs=example_uorfs)

    for uorf in calculations.uorfs:
        
        assert uorf in example_uorfs.five_prime_sequence

        assert uorf.startswith("ATG")
        endswith_stop_codon = uorf[-3:] in ["TAA", "TAG", "TGA"]
        assert endswith_stop_codon == True

        assert calculations.number_uorfs() == 3

        assert calculations.uorfs_lengths() == [51, 105, 66]

        assert calculations.gc_content() == [62.745098039215684, 70.47619047619048, 81.81818181818183]


def test_uorfs_with_20nt_more(example_uorfs: UORFs):

    longer_uorfs = UORFs_indel_analysis(uorfs= example_uorfs)
    
    result_20nt_diff = [
        (len(uorf_plus_20nt) - len(uorf)) <= 20 for uorf, uorf_plus_20nt in zip(example_uorfs.uorfs, longer_uorfs.uorfs_plus_20nt)
    ] 

    result_contains = [
        uorf in uorf_plus_20nt for uorf, uorf_plus_20nt in zip(example_uorfs.uorfs, longer_uorfs.uorfs_plus_20nt)
    ]

    result_length = [
        len(uorf) < len(uorf_plus_20nt) for uorf, uorf_plus_20nt in zip(example_uorfs.uorfs, longer_uorfs.uorfs_plus_20nt)
    ]

    assert result_20nt_diff == [True, True, True]
    assert result_contains == [True, True, True]
    assert result_length == [True, True, True]

def test_intercistonic_distances(example_uorfs: UORFs):
    
    assert UORFs_intercistonic_distance(uorfs= example_uorfs).intercistonic_distance_calculator() == [556, 216, 47]

def test_gc_content_10nt_after_uorf(example_uorfs: UORFs):

    assert UORFs_ten_nts_after(uorfs= example_uorfs).gc_content_10nt_after_uorf() == [90.0, 90.0, 70.0]
