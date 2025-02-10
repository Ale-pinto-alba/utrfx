import pytest

from utrfx.variant_util import prepare_alt_seq
from utrfx.model import FiveUTRCoordinates
from utrfx.genome import Contig, GenomicRegion, Strand, VariantCoordinates

class TestPrepareAltSeq:
    """
    A suite for testing `prepare_alt_seq` function.
    """

    CONTIG = Contig(
        name="X",
        gb_acc="GB_ACC",
        refseq_name="IRRELEVANT",
        ucsc_name="WHATEVER",
        length=100,
    )
    """
    A fake contig for testing
    """

    @pytest.fixture(scope="class")
    def five_utr_coordinates(self) -> FiveUTRCoordinates:
        return FiveUTRCoordinates(
            regions=(
                GenomicRegion(
                    TestPrepareAltSeq.CONTIG,
                    start=10,
                    end=25,
                    strand=Strand.POSITIVE,
                ),
                GenomicRegion(
                    TestPrepareAltSeq.CONTIG,
                    start=30,
                    end=35,
                    strand=Strand.POSITIVE,
                ),
            )
        )
    
    @pytest.fixture(scope="class")
    def pre_mrna_seq(self) -> str:
        """
        50 bases corresponding to a fake pre-mRNA sequence
        of a fake transcript.

        The sequence originates from the bases (10,50]
        of the `TestPrepareAltSeq.CONTIG`.
        """
        # Genomic coordinates (1-based):
        # 
        #      11       20        30        40        50
        #       |        |         |         |         |
        #       |                                      |
        #       |     5'UTR (1)    5'UTR (2)           |
        #       vvvvvvvvvvvvvvv     vvvvv              v
        return "AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"

    def test_snp(
        self,
        pre_mrna_seq: str,
        five_utr_coordinates: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "C", "T")

        #                    *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCTGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, pre_mrna_seq, five_utr_coordinates)
        assert actual == expected

    def test_del(
        self,
        pre_mrna_seq: str,
        five_utr_coordinates: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "CGG", "C")

        #                    ***
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, pre_mrna_seq, five_utr_coordinates)
        assert actual == expected

    def test_ins(
        self,
        pre_mrna_seq: str,
        five_utr_coordinates: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "C", "CTT")

        #                    *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCTTGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, pre_mrna_seq, five_utr_coordinates)
        assert actual == expected

    def test_mnv(
        self,
        pre_mrna_seq: str,
        five_utr_coordinates: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "CGG", "CA")

        #                    ***
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCAGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, pre_mrna_seq, five_utr_coordinates)
        assert actual == expected

    @staticmethod
    def make_variant(
        pos: int,
        ref: str,
        alt: str,
    ) -> VariantCoordinates:
        return VariantCoordinates.from_vcf_literal(
            TestPrepareAltSeq.CONTIG,
            pos=pos,
            ref=ref,
            alt=alt,
        )
