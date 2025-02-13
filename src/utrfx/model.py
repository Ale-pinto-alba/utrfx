import typing
import json
    
from utrfx.genome import GenomicRegion, Region


class FiveUTRCoordinates:
    """
    `FiveUTR` is a container for 5'UTR Genomic Regions.
    """
    def __init__(
        self,
        regions: typing.Collection[GenomicRegion],
    ):
        self._regions = tuple(regions)
        assert all(isinstance(r, GenomicRegion) for r in self._regions)
        # TODO: Low priority - check the regions are on the same strand

    @property
    def regions(self) -> typing.Collection[GenomicRegion]:
        return self._regions

    def __len__(self) -> int:
        return sum(len(region) for region in self._regions)
    
    def __repr__(self):
        regions_info = ", ".join([f"({region.contig.ucsc_name}, {region.start}, {region.end}, {region.strand})" for region in self._regions])
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
        assert isinstance(tx_id, str)
        self._tx_id = tx_id
        assert isinstance(five_utr, FiveUTRCoordinates)
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

    The UORF is upstream of the mORF.

    :param transcript: transcript with its corresponding identifier and 5'UTR Genomic Region(s).
    :param uorf: uORF region marked by its start and end nucleotide.
    :param ouorf: boolean indicating if the uORF overlaps with the mORF.
    """
    def __init__(
        self,
        five_utr: FiveUTRCoordinates,
        uorf: Region,
        ouorf: bool,
    ):
        assert isinstance(five_utr, FiveUTRCoordinates)
        self._five_utr = five_utr
        assert isinstance(uorf, Region)
        self._uorf = uorf
        assert isinstance(ouorf, bool)
        self._ouorf = ouorf

    @property
    def uorf(self) -> Region:
        return self._uorf
    
    @property
    def ouorf(self) -> bool:
        return self._ouorf

    def __len__(self) -> int:
        return len(self._uorf.end - self._uorf.start)  

    def __eq__(self, other):
        return (isinstance(other, UORFCoordinates)
                and self._five_utr== other._five_utr
                and self._uorf == other._uorf
                and self._ouorf == other._ouorf)
    
    def __repr__(self) -> str:
        return f"UORFCoordinates(Five_UTRs= {len(self._five_utr.regions)}, uORF= {self._uorf}, ouORF= {self._ouorf})"


class TxperGene:
    """
    `TxperGene` provides the canonical transcript of a given gene.
    """
    def __init__(
        self,
        fpath: str,
    ):
        self._fpath = fpath
        self._dict = self.obtain_tx_as_dictionary_from_json(self._fpath)

    @staticmethod
    def obtain_tx_as_dictionary_from_json(
        fpath: str,
    ) -> dict:
        """
        Load transcript data from a JSON file into a dictionary.
        """
        try:
            with open(fpath, "r", encoding="utf-8") as file:
                return json.load(file)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            raise ValueError(f"Error: {e}")

    def ensembl_transcript(
        self,
        gene_symbol: str,
    ) -> typing.Optional[str]:
        """
        Retrieve the canonical ENSEMBL transcript if available.
        """
        return self._dict.get(gene_symbol)