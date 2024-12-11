import pytest
import typing

from utrfx.genome import GenomeBuild, GenomicRegion, Region, GRCh38, Strand
from utrfx.model import FiveUTR, UORF
from utrfx.uorf import download_fasta_from_ensembl, get_five_prime_sequence, uorf_extractor, gc_content, gc_content_n_bases_downstream, uorfs_plus_n_nts_downstream_extractor, intercistonic_distance

@pytest.fixture(scope="module")
def transcript_fasta() -> str:
    return download_fasta_from_ensembl("ENST00000381418")


@pytest.fixture(scope="module")
def five_utr() -> FiveUTR:
    return FiveUTR([GenomicRegion(contig= GenomeBuild.contig_by_name(GRCh38, "8"), start= 22114419, end= 22115043, strand= Strand.NEGATIVE)])


@pytest.fixture(scope="module")
def five_utr_sequence(transcript_fasta: str, five_utr: FiveUTR) -> str:
    return get_five_prime_sequence(transcript_sequence= transcript_fasta, five_utrs= five_utr)


def test_uorf_extractor(five_utr: FiveUTR, five_utr_sequence: str) -> typing.Collection[UORF]:
    uorfs = uorf_extractor(five_utrs= five_utr, five_utr_seq= five_utr_sequence)
    assert len(uorfs) == 3

    first_uorf, second_uorf, third_uorf = uorfs

    assert first_uorf._uorf.start == 16
    assert first_uorf._uorf.end == 67

    assert second_uorf._uorf.start == 302
    assert second_uorf._uorf.end == 407

    assert third_uorf._uorf.start == 510
    assert third_uorf._uorf.end == 576


@pytest.mark.parametrize(
        "uorf, expected",
        [
            (UORF(five_utr= five_utr, uorf= Region(start= 16, end= 67)), 62.745098039215684),
            (UORF(five_utr= five_utr, uorf= Region(start= 302, end= 407)), 70.47619047619048),
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 81.81818181818183),
        ]
)
def test_gc_content(five_utr_sequence: str, uorf, expected):
    assert gc_content(five_utr_seq= five_utr_sequence, uorf= uorf) == expected


def test_uorf_ends_out_of_five_prime(five_utr_sequence: str):
    overlapping_uorf = UORF(five_utr= five_utr, uorf= Region(start= 700, end= 800))

    with pytest.raises(ValueError) as e:
        gc_content(five_utr_seq= five_utr_sequence, uorf= overlapping_uorf)

    assert e.value.args == ("uORF overlaps with the mORF",)


@pytest.mark.parametrize(
        "uorf, bases, expected",
        [
            (UORF(five_utr= five_utr, uorf= Region(start= 16, end= 67)), 10, 90.0),
            (UORF(five_utr= five_utr, uorf= Region(start= 302, end= 407)), 10, 90.0),
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 10, 60.0),
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 700, 72.91666666666666), # Number of bases out of 5'UTR region,
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 8500, 72.91666666666666), # so the GC content remains the same, clipped to 48 bases
        ]
)
def test_gc_content_n_bases_downstream(five_utr_sequence: str, uorf, bases, expected):
    assert gc_content_n_bases_downstream(five_utr_seq= five_utr_sequence, uorf= uorf, bases= bases) == expected


@pytest.mark.parametrize(
        "uorf, bases, expected",
        [
            (UORF(five_utr= five_utr, uorf= Region(start= 16, end= 67)), 10, 61),
            (UORF(five_utr= five_utr, uorf= Region(start= 302, end= 407)), 10, 115),
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 10, 76),
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 700, 114), 
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 8500, 114),
        ]
)
def test_uorfs_plus_nts_downstream(five_utr_sequence: str, uorf, bases, expected):
    assert len(uorfs_plus_n_nts_downstream_extractor(five_utr_seq= five_utr_sequence, uorf= uorf, bases= bases)) == expected


@pytest.mark.parametrize(
        "uorf, expected",
        [
            (UORF(five_utr= five_utr, uorf= Region(start= 16, end= 67)), 557),
            (UORF(five_utr= five_utr, uorf= Region(start= 302, end= 407)), 217),
            (UORF(five_utr= five_utr, uorf= Region(start= 510, end= 576)), 48),
        ]
)
def test_intercistonic_distances(five_utr_sequence: str, uorf, expected):
    assert intercistonic_distance(five_utr_seq= five_utr_sequence, uorf= uorf) == expected