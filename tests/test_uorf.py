import pytest
import typing
import os

from utrfx.uorf import UORFsProcessor

def test_utr_length(example_uorfs: UORFsProcessor):
    assert example_uorfs.five_utr_lenght() == 623

def test_uorfs(example_uorfs: UORFsProcessor):
    for uorf in example_uorfs._uorfs:
        assert uorf in example_uorfs.tx_sequence         # Check each uORF is in the transcript sequence               
        assert uorf in example_uorfs.five_utr_sequence   # Check each uORF is in the 5'UTR sequence
        assert uorf.startswith("ATG")
        endswith_stop_codon = uorf[-3:] in ["TAA", "TAG", "TGA"]
        assert endswith_stop_codon == True

def test_uorfs_with_20nt_more(example_uorfs: UORFsProcessor):
    
    result_20nt_diff = [
        (len(uorf_plus_20nt) - len(uorf)) <= 20 for uorf, uorf_plus_20nt in zip(example_uorfs.uorfs, example_uorfs.uorfs_with_20nt_more)
    ] 

    result_contains = [
        uorf in uorf_plus_20nt for uorf, uorf_plus_20nt in zip(example_uorfs.uorfs, example_uorfs.uorfs_with_20nt_more)
    ]

    result_length = [
        len(uorf) < len(uorf_plus_20nt) for uorf, uorf_plus_20nt in zip(example_uorfs.uorfs, example_uorfs.uorfs_with_20nt_more)
    ]

    assert result_20nt_diff == [True, True, True]
    assert result_contains == [True, True, True]
    assert result_length == [True, True, True]

def test_number_uorfs(example_uorfs: UORFsProcessor):
    
    assert example_uorfs.tx_id == "ENST00000381418.9"
    assert example_uorfs.number_of_uorfs() == 3

def test_uorfs_lengths(example_uorfs: UORFsProcessor):
    
    assert example_uorfs.uorfs_lengths() == [51, 105, 66]

def test_gc_content(example_uorfs: UORFsProcessor):

    assert example_uorfs.gc_content() == [62.745098039215684, 70.47619047619048, 81.81818181818183]

def test_intercistonic_distances(example_uorfs: UORFsProcessor):
    
    assert example_uorfs.intercistonic_distance() == [556, 216, 47]

def test_gc_content_10nt_after_uorf(example_uorfs: UORFsProcessor):

    assert example_uorfs.gc_content_10nt_after_uorf() == [90.0, 90.0, 70.0]

