import pytest

from utrfx.genome import Region
from utrfx.model import FiveUTRCoordinates, UORFCoordinates
from utrfx.uorf import gc_content, gc_content_n_bases_downstream, uorfs_plus_n_nts_downstream_extractor, intercistronic_distance


@pytest.mark.parametrize(
    "region, expected",
    [
        (Region(start=16, end=67), ((13+19)/51)),
        (Region(start=302, end=407), ((31+43)/105)),
        (Region(start=510, end=576), ((26+28)/66)),
    ]
)
def test_gc_content(
    hbb_five_utr_sequence: str, 
    hbb_five_utr: FiveUTRCoordinates, 
    region: Region, 
    expected: float,
):
    uorf = UORFCoordinates(five_utr=hbb_five_utr, uorf=region)

    assert gc_content(five_sequence=hbb_five_utr_sequence, uorf=uorf) == expected


def test_uorf_ends_out_of_five_prime(
    hbb_five_utr_sequence: str,
    hbb_five_utr: FiveUTRCoordinates,
):
    overlapping_uorf = UORFCoordinates(five_utr=hbb_five_utr, uorf=Region(start=700, end=800))

    with pytest.raises(ValueError) as e:
        gc_content(five_sequence=hbb_five_utr_sequence, uorf=overlapping_uorf)

    assert e.value.args == ("uORF overlaps with the mORF",)


@pytest.mark.parametrize(
        "region, bases, expected",
        [
            ((Region(start=16, end=67)), 10, ((4+5)/10)),
            ((Region(start=302, end=407)), 10, ((3+6)/10)),
            ((Region(start=510, end=576)), 10, ((4+3)/10)),
            ((Region(start=510, end=576)), 700, ((18+15)/47)), # Number of bases out of 5'UTR region,
            ((Region(start=510, end=576)), 8500, ((18+15)/47)), # so the GC content remains the same, clipped to 48 bases
        ]
)
def test_gc_content_n_bases_downstream(
    hbb_five_utr_sequence: str,
    hbb_five_utr: FiveUTRCoordinates,
    region: Region,
    bases: int,
    expected: float, 
):
    uorf = UORFCoordinates(five_utr=hbb_five_utr, uorf=region)

    assert gc_content_n_bases_downstream(five_sequence=hbb_five_utr_sequence, uorf=uorf, bases=bases) == expected


@pytest.mark.parametrize(
        "region, bases, expected",
        [
            ((Region(start=16, end=67)), 10, 61),
            ((Region(start=302, end=407)), 10, 115),
            ((Region(start=510, end=576)), 10, 76),
            ((Region(start=510, end=576)), 700, 113), 
            ((Region(start=510, end=576)), 8500, 113),
        ]
)
def test_uorfs_plus_nts_downstream(
    hbb_five_utr_sequence: str,
    hbb_five_utr: FiveUTRCoordinates,
    region: Region,
    bases: int,
    expected: float,
):
    uorf = UORFCoordinates(five_utr=hbb_five_utr, uorf=region)

    assert len(uorfs_plus_n_nts_downstream_extractor(five_sequence=hbb_five_utr_sequence, uorf=uorf, bases=bases)) == expected


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), 556),
            ((Region(start=302, end=407)), 216),
            ((Region(start=510, end=576)), 47),
        ]
)
def test_intercistonic_distances(
    hbb_five_utr_sequence: str,
    hbb_five_utr: FiveUTRCoordinates,
    region: Region,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hbb_five_utr, uorf=region)

    assert intercistronic_distance(five_sequence=hbb_five_utr_sequence, uorf=uorf) == expected