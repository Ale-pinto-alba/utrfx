import pytest

from utrfx.genome import Region
from utrfx.model import FiveUTRCoordinates, UORFCoordinates
from utrfx.uorf import gc_content, gc_content_n_bases_downstream, uorfs_plus_n_nts_downstream_extractor, intercistronic_distance, cap_five_to_uorf_distance, kozak_sequence_strength


@pytest.mark.parametrize(
    "region, expected",
    [
        (Region(start=16, end=67), ((13+19)/51)),
        (Region(start=302, end=407), ((31+43)/105)),
        (Region(start=510, end=576), ((26+28)/66)),
    ]
)
def test_gc_content(
    hr_five_utr_sequence: str, 
    hr_five_utr: FiveUTRCoordinates, 
    region: Region, 
    expected: float,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf= False)
    result = gc_content(five_sequence=hr_five_utr_sequence, uorf=uorf)
    
    assert result == pytest.approx(expected, rel=1e-6)


def test_uorf_ends_out_of_five_prime(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
):
    overlapping_uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=Region(start=700, end=800), ouorf=True)

    with pytest.raises(ValueError) as e:
        gc_content(five_sequence=hr_five_utr_sequence, uorf=overlapping_uorf)

    assert e.value.args == ("uORF overlaps with the mORF",)


@pytest.mark.parametrize(
        "region, bases, ouorf, expected",
        [
            ((Region(start=16, end=67)), 10, False, ((4+5)/10)),
            ((Region(start=302, end=407)), 10, False, ((3+6)/10)),
            ((Region(start=510, end=576)), 10, False, ((4+3)/10)),
            ((Region(start=510, end=576)), 700, True, ((18+15)/47)), # Number of bases out of 5'UTR region,
            ((Region(start=510, end=576)), 8500, True, ((18+15)/47)), # so the GC content remains the same, clipped to 47 bases
        ]
)
def test_gc_content_n_bases_downstream(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    bases: int,
    ouorf: bool,
    expected: float, 
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)
    result = gc_content_n_bases_downstream(five_sequence=hr_five_utr_sequence, uorf=uorf, bases=bases) 

    assert result == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize(
        "region, bases, ouorf, expected",
        [
            ((Region(start=16, end=67)), 10, False, 61),
            ((Region(start=302, end=407)), 10, False, 115),
            ((Region(start=510, end=576)), 10, False, 76),
            ((Region(start=510, end=576)), 700, True, 113), 
            ((Region(start=510, end=576)), 8500, True, 113),
        ]
)
def test_uorfs_plus_nts_downstream(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    bases: int,
    ouorf: bool,
    expected: float,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)

    assert len(uorfs_plus_n_nts_downstream_extractor(five_sequence=hr_five_utr_sequence, uorf=uorf, bases=bases)) == expected


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), 556),
            ((Region(start=302, end=407)), 216),
            ((Region(start=510, end=576)), 47),
        ]
)
def test_intercistronic_distances(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=False)

    assert intercistronic_distance(five_sequence=hr_five_utr_sequence, uorf=uorf) == expected


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), 16),
            ((Region(start=302, end=407)), 302),
            ((Region(start=510, end=576)), 510),
        ]
)
def test_five_cap_to_uorf_distance(
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=False)

    assert cap_five_to_uorf_distance(uorf=uorf) == expected


@pytest.mark.parametrize(
        "region, ouorf, expected",
        [
            ((Region(start=16, end=67)), False, 2),
            ((Region(start=302, end=407)), False, 1),
            ((Region(start=510, end=576)), False, 1),
            ((Region(start=606, end=623)), True, 0)
        ]
)
def test_kozak_sequence_strength(
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
    region: Region,
    ouorf: bool,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hr_five_utr, uorf=region, ouorf=ouorf)

    assert kozak_sequence_strength(five_sequence=hr_five_utr_sequence, uorf=uorf) == expected