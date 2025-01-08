import pytest

from utrfx.genome import Region
from utrfx.model import FiveUTRCoordinates, UORFCoordinates
from utrfx.uorf import gc_content, gc_content_n_bases_downstream, uorfs_plus_n_nts_downstream_extractor, intercistronic_distance, cap_five_to_uorf_distance


@pytest.mark.parametrize(
    "region, expected",
    [
        (Region(start=16, end=67), 0.62745098039215684),
        (Region(start=302, end=407), 0.7047619047619048),
        (Region(start=510, end=576), 0.8181818181818182),
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
            ((Region(start=16, end=67)), 10, 0.9),
            ((Region(start=302, end=407)), 10, 0.9),
            ((Region(start=510, end=576)), 10, 0.6),
            ((Region(start=510, end=576)), 700, 0.7291666666666666), # Number of bases out of 5'UTR region,
            ((Region(start=510, end=576)), 8500, 0.7291666666666666), # so the GC content remains the same, clipped to 48 bases
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
            ((Region(start=510, end=576)), 700, 114), 
            ((Region(start=510, end=576)), 8500, 114),
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
            ((Region(start=16, end=67)), 557),
            ((Region(start=302, end=407)), 217),
            ((Region(start=510, end=576)), 48),
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


@pytest.mark.parametrize(
        "region, expected",
        [
            ((Region(start=16, end=67)), 16),
            ((Region(start=302, end=407)), 302),
            ((Region(start=510, end=576)), 510),
        ]
)
def test_five_cap_to_uorf_distance(
    hbb_five_utr: FiveUTRCoordinates,
    region: Region,
    expected: int,
):
    uorf = UORFCoordinates(five_utr=hbb_five_utr, uorf=region)

    assert cap_five_to_uorf_distance(uorf=uorf) == expected