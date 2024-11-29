import typing
    
from .genome import GenomicRegion

class FiveUTR:
    """
    `FiveUTR` is a container for 5'UTR Genomic Regions.
    """

    def __init__(
        self,
        regions: typing.Collection[GenomicRegion],
    ):
        self._regions = regions

    def __len__(self) -> int:
        return 0 if len(self._regions) == 0 else sum(len(region) for region in self._regions)
    
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
    