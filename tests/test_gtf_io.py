import os

import pytest

from utrfx.genome import GenomeBuild, Strand
from utrfx.gtf_io import read_gtf_into_txs

class TestGtfIo:

    @pytest.fixture(scope="class")
    def fpath_example_gtf(self, fpath_data_dir: str) -> str:  
        return os.path.join(fpath_data_dir,  "Homo.sapiens.GRCh38_sample.gtf")

    def test_read_gtf_into_txs(
        self,
        fpath_example_gtf: str,
        genome_build: GenomeBuild,
    ):
        transcripts = read_gtf_into_txs(fpath_example_gtf, genome_build)

        assert len(transcripts) == 1_327

        # Positive strand
        for tx in transcripts:
            if tx.tx_id == "ENST00000432186.6":
                our_favorite_tx = tx
                break
        
        assert our_favorite_tx is not None
        assert len(our_favorite_tx.five_utr.regions) == 2
        one, two = sorted(our_favorite_tx.five_utr.regions, key=lambda region: region.start)

        assert one.contig.name == "22"
        assert one.start == 44668712
        assert one.end == 44668805
        assert one.strand == Strand.POSITIVE
        assert two.contig.name == "22"
        assert two.start == 44702491
        assert two.end == 44702501
        assert two.strand == Strand.POSITIVE

        # Negative strand
        for tx in transcripts:
            if tx.tx_id == "ENST00000703965.1":
                our_another_favorite_tx = tx
                break

        assert our_another_favorite_tx is not None
        assert len(our_another_favorite_tx._five_utr.regions) == 2
        three, four = sorted(our_another_favorite_tx.five_utr.regions, key=lambda region: region.start)

        assert three.contig.name == "22"
        assert three.start == 26837999
        assert three.end == 26838057
        assert three.strand == Strand.NEGATIVE
        assert four.contig.name == "22"
        assert four.start == 26841401
        assert four.end == 26841576
        assert four.strand == Strand.NEGATIVE