import typing
    
from .genome import GenomicRegion

class Region:
    """
    
    Hey, please use static constructors to create a region.

    >>> from utrfx.genome import Region
    >>> a = Region.from_one_based("chr1", 1, 10)
    >>> a.length()
    10
    """

    @staticmethod
    def from_one_based(contig: str, start: int, end: int) -> "Region":
        """
        :param start: one based start coordinate (inclusive)
        """
        return Region(contig=contig, start=start-1, end=end)
    
    def __init__(
        self,
        contig: str,
        start: int,
        end: int,
    ):
        assert isinstance(contig, str), "contig given is not a str"
        self._contig = contig

        assert isinstance(start, int), "start is not an int"
        self._start = start
        assert self._start >= 0, "start must be non-negative"

        assert isinstance(end, int), "end is not an int"
        self._end = end
        assert self._end >= 0, "end must be non-negative"

        assert self._start <= self._end, "end cannot be before start"

    def contig(self) -> str:
        return self._contig
    
    def length(self) -> int:
        return self._end - self._start

    def overlaps_with(self, other: "Region") -> bool:

        if self._contig == other._contig: 
            if other._start == self._start and other._end == self._end:
                return True
            elif other._start <= self._end and other._end > self._end:
                return True
            elif other._end >= self._start and other._start < self._start:
                return True
            else:
                return False
        else:
            return False


class FiveUTR:
    """
    `FiveUTR` is a container for 5'UTR Genomic Regions.
    """

    def __init__(
        self,
        regions: typing.Collection[GenomicRegion],
    ):
        self._regions = regions
    
    def __repr__(self):
        regions_info = ", ".join([f"({region._contig.ucsc_name}, {region.start}, {region.end}, {region.strand})" for region in self._regions])
        return f"FiveUTR(regions={len(self._regions)} regions: {regions_info})"

    @property
    def regions(self) -> typing.Collection[GenomicRegion]:
        return self._regions
    

class Transcript:
    """
    `Transcript` represents the 5'UTR Genomic Region(s) of a transcript.
    """
    
    def __init__(
        self, 
        tx_id: str, 
        five_utr: FiveUTR
    ):
        self._tx_id = tx_id
        self._five_utr = five_utr

    @property
    def tx_id(self) -> str:
        return self._tx_id
    
    @property
    def five_utr(self) -> FiveUTR:
        return self._five_utr

    def __repr__(self):
        return f"Transcript(tx_id={self._tx_id}, five_utr={repr(self._five_utr)})"