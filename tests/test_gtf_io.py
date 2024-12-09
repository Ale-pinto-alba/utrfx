import os

import pytest

from utrfx.genome import GenomeBuild, Strand, Region, GenomicRegion, Contig
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
        assert one.start == 44_668_712
        assert one.end == 44_668_805
        assert one.strand == Strand.POSITIVE
        assert two.contig.name == "22"
        assert two.start == 44_702_491
        assert two.end == 44_702_501
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
        assert three.start == 26_837_999
        assert three.end == 26_838_057
        assert three.strand == Strand.NEGATIVE
        assert four.contig.name == "22"
        assert four.start == 26_841_401
        assert four.end == 26_841_576
        assert four.strand == Strand.NEGATIVE

        three_in_positive_strand = three.with_strand(other= Strand.POSITIVE)
        three_one_based_start = three_in_positive_strand.start + 1
        assert three_one_based_start == 23_980_412
        assert three_in_positive_strand.end == 23_980_469

        four_in_positive_strand = three.with_strand(other= Strand.POSITIVE)
        four_one_based_start = four_in_positive_strand.start + 1
        assert four_one_based_start == 23_980_412
        assert four_in_positive_strand.end == 23_980_469
        