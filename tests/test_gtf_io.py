import pytest

from utrfx.genome import GenomeBuild, GenomicRegion, Contig, Strand
from utrfx.gtf_io import read_gtf_into_txs

class TestGtfIo:

    def test_read_gtf_into_txs(
        self,
        fpath_example_gtf: str,
        genome_build: GenomeBuild,
        first_region_favorite_tx,
        second_region_favorite_tx,
        first_region_another_tx,
        second_region_another_tx,
    ):
        transcripts = read_gtf_into_txs(fpath_example_gtf, genome_build)

        assert len(transcripts) == 1_327

        # Positive strand
        for tx in transcripts:
            if tx.tx_id == "ENST00000432186.6":
                our_favorite_tx = tx
                break
        
        assert our_favorite_tx is not None
        assert len(our_favorite_tx._five_utr._regions) == 2 
        assert our_favorite_tx._five_utr._regions[0] == first_region_favorite_tx
        assert our_favorite_tx._five_utr._regions[1] == second_region_favorite_tx

        # Negative strand
        for tx in transcripts:
            if tx.tx_id == "ENST00000703965.1":
                our_another_favorite_tx = tx
                break

        assert our_another_favorite_tx is not None
        assert len(our_another_favorite_tx._five_utr._regions) == 2
        assert our_another_favorite_tx._five_utr._regions[0] == first_region_another_tx
        assert our_another_favorite_tx._five_utr._regions[1] == second_region_another_tx