import abc

from utrfx.genome import GenomicRegion


class GenomicSequenceService(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def fetch(self, query: GenomicRegion) -> str:
        pass


class PySamService(GenomicSequenceService):

    def __init__(
        self,
        fpath_fasta: str,
    ):
        # TODO - use PySAM to open an indexed FASTA file for random sequence retrieval
        self._fasta_handle = None  # replace with an actual fasta file
        pass

    def fetch(self, query: GenomicRegion) -> str:
        # get the coordinates on the Positive strand 
        # fetch sequence
        # reverse complement sequence if `query.strand == -`
        raise NotImplementedError
