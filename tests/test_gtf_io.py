import pytest

from utrfx.genome import GenomeBuild, GenomicRegion, Contig, Strand
from utrfx.gtf_io import read_gtf_into_txs

class TestGtfIo:
    
    @pytest.fixture(scope="class")
    def first_region_favorite_tx(self, genome_build: GenomeBuild):
        contig = genome_build.contig_by_name("chr22")
        return GenomicRegion(contig=contig, start=44668712, end=44668805, strand=Strand.POSITIVE)

    @pytest.fixture(scope="class")
    def second_region_favorite_tx(self, genome_build: GenomeBuild):
        contig = genome_build.contig_by_name("chr22")
        return GenomicRegion(contig=contig, start=44702491, end=44702501, strand=Strand.POSITIVE)

    @pytest.fixture(scope="class")
    def first_region_another_tx(self, genome_build: GenomeBuild):
        contig = genome_build.contig_by_name("chr22")
        return GenomicRegion(contig=contig, start=26837999, end=26838057, strand=Strand.NEGATIVE)

    @pytest.fixture(scope="class")
    def second_region_another_tx(self, genome_build: GenomeBuild):
        contig = genome_build.contig_by_name("chr22")
        return GenomicRegion(contig=contig, start=26841401, end=26841576, strand=Strand.NEGATIVE)

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
        assert len(our_favorite_tx._five_utr.regions) == 2 
        assert our_favorite_tx._five_utr.regions[0].contig.name == first_region_favorite_tx.contig.name
        assert our_favorite_tx._five_utr.regions[0].start == first_region_favorite_tx.start
        assert our_favorite_tx._five_utr.regions[0].end == first_region_favorite_tx.end
        assert our_favorite_tx._five_utr.regions[0].strand == first_region_favorite_tx.strand
        assert our_favorite_tx._five_utr.regions[1].contig.name == second_region_favorite_tx.contig.name
        assert our_favorite_tx._five_utr.regions[1].start == second_region_favorite_tx.start
        assert our_favorite_tx._five_utr.regions[1].end == second_region_favorite_tx.end
        assert our_favorite_tx._five_utr.regions[1].strand == second_region_favorite_tx.strand

        # Negative strand
        for tx in transcripts:
            if tx.tx_id == "ENST00000703965.1":
                our_another_favorite_tx = tx
                break

        assert our_another_favorite_tx is not None
        assert len(our_another_favorite_tx._five_utr.regions) == 2
        assert our_another_favorite_tx._five_utr.regions[0].contig.name == first_region_another_tx.contig.name
        assert our_another_favorite_tx._five_utr.regions[0].start == first_region_another_tx.start
        assert our_another_favorite_tx._five_utr.regions[0].end == first_region_another_tx.end
        assert our_another_favorite_tx._five_utr.regions[0].strand == first_region_another_tx.strand
        assert our_another_favorite_tx._five_utr.regions[1].contig.name == second_region_another_tx.contig.name
        assert our_another_favorite_tx._five_utr.regions[1].start == second_region_another_tx.start
        assert our_another_favorite_tx._five_utr.regions[1].end == second_region_another_tx.end
        assert our_another_favorite_tx._five_utr.regions[1].strand == second_region_another_tx.strand