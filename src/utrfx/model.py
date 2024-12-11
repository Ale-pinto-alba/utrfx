import typing
    
from utrfx.genome import GenomicRegion, Region


class FiveUTRCoordinates:
    """
    `FiveUTRCoordinates` is a container for 5'UTR Genomic Regions.

    The length of FiveUTR corresponds to the sum of lengths of all its regions.

    We "kind of" assume all regions are on the same strand.

    :param regions: list of Genomic Regions corresponding to the 5'UTRs regions.
    """
    def __init__(
        self,
        regions: typing.Collection[GenomicRegion],
    ):
        # TODO: Low priority - check the regions are on the same strand
        self._regions = regions

    @property
    def regions(self) -> typing.Collection[GenomicRegion]:
        return self._regions

    def __len__(self) -> int:
        return sum(len(region) for region in self._regions)
    
    def __repr__(self):
        regions_info = ", ".join([f"({region._contig.ucsc_name}, {region.start}, {region.end}, {region.strand})" for region in self._regions])
        return f"FiveUTRCoordinates(regions={len(self._regions)} regions: {regions_info})"
    

class TranscriptCoordinates:
    """
    `TranscriptCoordinates` represents the 5'UTR Genomic Region(s) of a transcript.

    :param tx_id: transcript identifier, e.g. `ENST00000381418`
    :param five_utr: 5'UTR Genomic Region(s). 
    """
    def __init__(
        self, 
        tx_id: str, 
        five_utr: FiveUTRCoordinates
    ):
        self._tx_id = tx_id
        self._five_utr = five_utr

    @property
    def tx_id(self) -> str:
        return self._tx_id
    
    @property
    def five_utr(self) -> FiveUTRCoordinates:
        return self._five_utr

    def __repr__(self):
        return f"TranscriptCoordinates(tx_id={self._tx_id}, five_utr={repr(self._five_utr)})"


class UORFCoordinates:
    """
    `UORFCoordinates` represents an uORF of a transcript.

    The UORF is upstream of the mORF and they do *not* overlap.

    :param transcript: transcript with its corresponding identifier and 5'UTR Genomic Region(s).
    :param uorf: uORF region marked by its start and end nucleotide.
    """
    def __init__(
        self,
        five_utr: FiveUTRCoordinates,
        uorf: Region,
    ):
        self._five_utr = five_utr
        self._uorf = uorf

    @property
    def uorf(self) -> Region:
        return self._uorf

    def __len__(self) -> int:
        return len(self._uorf.end - self._uorf.start)  

    def __eq__(self, other):
        return (isinstance(other, UORFCoordinates)
                and self._five_utr== other._five_utr
                and self._uorf == other._uorf)
    
    def __repr__(self) -> str:
        return f"UORFCoordinates(Five_UTRs= {len(self._five_utr.regions)}, uORF= {self._uorf})"