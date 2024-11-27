import pytest

from utrfx.model import Region

@pytest.mark.parametrize(
    "start, end, expected",
    [
        (-5, 10, "start must be non-negative"),
        (5, -10, "end must be non-negative"),
        (10, 9, "end cannot be before start"),
    ]
)
def test_explodes_with_invalid_coordinates(
    start: int,
    end: int,
    expected: str,
):
    with pytest.raises(AssertionError) as e:
        Region(contig="whatever", start=start, end=end)

    assert e.value.args == (expected,)

@pytest.mark.parametrize(
    "start, end, expected",
    [
        (3.4, 8, "start is not an int"),
        (3+4j, 5, "start is not an int"),
        ("No_vale", 4, "start is not an int"),
        (8, 8+5j, "end is not an int"),
        (8, "No_vale", "end is not an int"),
        (8, 10.4, "end is not an int"),    
    ]
)
def test_invalid_chromosomes(
    start: int,
    end: int,
    expected: str,
):
    with pytest.raises(AssertionError) as e:
        Region(contig="chr1", start=start, end=end)

    assert e.value.args == (expected,)

@pytest.mark.parametrize(
        "contig, expected",
        [
            ("chr1", "chr1"),
            ("chr8", "chr8"),
            ("chrX", "chrX"),
        ]
)
def test_contig(contig, expected):
    assert Region(contig=contig, start=1, end=5).contig() == expected

@pytest.mark.parametrize(
    "start, end, expected",
    [
        (1, 10, 10),
        (15, 115, 101),
        (8, 9, 2)
    ]
)
def test_length(start, end, expected):
    assert Region.from_one_based(contig="chr1", start=start, end=end).length() == expected

@pytest.mark.parametrize(
    "region, other, expected",                        
    [
        (Region.from_one_based("chr2", 10, 20), Region.from_one_based("chr1", 1, 10), False),
        (Region.from_one_based("chr1", 10, 20), Region.from_one_based("chr1", 115, 205), False),
        (Region.from_one_based("chr1", 100, 250), Region.from_one_based("chr1", 4, 5), False),
        (Region.from_one_based("chr1", 1, 10), Region.from_one_based("chr1", 1, 10), True),
        (Region.from_one_based("chr2", 10, 20), Region.from_one_based("chr2", 15, 30), True),
        (Region.from_one_based("chr1", 10, 20), Region.from_one_based("chr1", 19, 30), True),
    ]
)
def test_overlaps_with(region, other, expected):
    assert region.overlaps_with(other) == expected
