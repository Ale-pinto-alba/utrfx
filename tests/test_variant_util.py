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
    def five_utr_coordinates_forward(self) -> FiveUTRCoordinates:
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
                    end=55,
                    strand=Strand.POSITIVE,
                ),
            )
        )
    
    @pytest.fixture(scope="class")
    def cdna_seq(self) -> str:
        """
        40 bases corresponding to a fake cDNA sequence
        of the 5'UTR region of a fake transcript.

        The sequence originates from the bases (10,50]
        of the `TestPrepareAltSeq.CONTIG`.
        """
        # Genomic coordinates (1-based):
        # 
        #      11            25    31            45        55
        #       |             |     |             |         |
        #       |             |     |                       |
        #       |  5'UTR (1)  |     |       5'UTR (2)       |
        #       v             v     v                       v
        return "AAAAACCCCCGGGGG" + "TTTTTAAAAACCCCCGGGGGTTTTT"

    def test_snp(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "C", "T")

        #                    *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCTGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_del(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "CGG", "C")

        #                    ***
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_ins(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "C", "CTT")

        #                    *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCTTGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_mnv(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(20, "CGG", "CA")

        #                    ***
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCAGGGTTTTTAAAAACCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_snp_falls_on_second_five_utr_region(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(40, "A", "T")

        #                                   *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCGGGGGTTTTTAAAATCCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_del_falls_on_second_five_utr_region(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(40, "ACC", "A")

        #                                   *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCGGGGGTTTTTAAAAACCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_ins_falls_on_second_five_utr_region(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(40, "A", "ATT")

        #                                   *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCGGGGGTTTTTAAAAATTCCCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    def test_mnv_falls_on_second_five_utr_region(
        self,
        cdna_seq: str,
        five_utr_coordinates_forward: FiveUTRCoordinates,
    ):
        vc = TestPrepareAltSeq.make_variant(40, "ACC", "AT")

        #                                   *
        #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
        expected = "AAAAACCCCCGGGGGTTTTTAAAAATCCCGGGGGTTTTT"
        actual = prepare_alt_seq(vc, cdna_seq, five_utr_coordinates_forward)
        assert actual == expected

    class TestNegativeStrand:

        @pytest.fixture(scope="class")
        def five_utr_coordinates_negative(self) -> FiveUTRCoordinates:
            return FiveUTRCoordinates(
                regions=(
                    GenomicRegion(
                        TestPrepareAltSeq.CONTIG,
                        start=5,
                        end=15,
                        strand=Strand.NEGATIVE,
                    ),
                    GenomicRegion(
                        TestPrepareAltSeq.CONTIG,
                        start=20,
                        end=35,
                        strand=Strand.NEGATIVE,
                    ),
                )
            )

        @pytest.fixture(scope="class")
        def cdna_seq_negative(self) -> str:
            """
            25 bases corresponding to a fake cDNA sequence
            of the 5'UTR region of a fake transcript on the negative strand.
            The sequence originates from the bases
            spanned by (5,15](-) (20,35](-)
            regions of the `TestPrepareAltSeq.CONTIG`.
            """
            # Genomic coordinates (1-based):
            # 
            #        5'UTR (1)        5'UTR (2)
            #       6       15    21            35
            #       |        |     |             |
            #       |        |     |             |
            #       |        |     |             |
            #       v        v     v             v
            return "GGGGGCCCCC" + "TTTTTAAAAACCCCC"

        def test_snp(
            self,
            cdna_seq_negative: str,
            five_utr_coordinates_negative: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(91, "C", "A")

            #               *
            #     ref:  GGGGGCCCCCTTTTTAAAAACCCCC
            expected = "GGGGTCCCCCTTTTTAAAAACCCCC"
            actual = prepare_alt_seq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            assert actual == expected

        def test_del(
            self,
            cdna_seq_negative: str,
            five_utr_coordinates_negative: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(91, "CC", "C")

            #              **
            #     ref:  GGGGGCCCCCTTTTTAAAAACCCCC
            expected = "GGGGCCCCCTTTTTAAAAACCCCC"
            actual = prepare_alt_seq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            assert actual == expected

        def test_ins(
            self,
            cdna_seq_negative: str,
            five_utr_coordinates_negative: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(91, "C", "CA")

            #               *
            #     ref:  GGGGGCCCCCTTTTTAAAAACCCCC
            expected = "GGGGGTCCCCCTTTTTAAAAACCCCC"
            actual = prepare_alt_seq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            assert actual == expected

        def test_mnv(
            self,
            cdna_seq_negative: str,
            five_utr_coordinates_negative: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(70, "TG", "GAC")

            #                              **
            #     ref:  GGGGGCCCCCTTTTTAAAAACCCCC
            expected = "GGGGGCCCCCTTTTTAAAACTGCCCC"
            actual = prepare_alt_seq(vc, cdna_seq_negative, five_utr_coordinates_negative)
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
    
    
def test_hr_variant_one(
    hr_variant_one: VariantCoordinates,
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,
):
    vc = prepare_alt_seq(hr_variant_one, hr_five_utr_sequence, hr_five_utr)

    assert hr_five_utr_sequence[405] == "A"
    assert vc[405] == "G"

def test_hr_variant_two(
    hr_variant_two: VariantCoordinates,
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,        
):
    vc = prepare_alt_seq(hr_variant_two, hr_five_utr_sequence, hr_five_utr)

    assert hr_five_utr_sequence[321] == "C"
    assert vc[321] == "A"

def test_hr_variant_three(
    hr_variant_three: VariantCoordinates,
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,        
):
    vc = prepare_alt_seq(hr_variant_three, hr_five_utr_sequence, hr_five_utr)

    assert hr_five_utr_sequence[308] == "C"
    assert vc[308] == "T"

def test_hr_variant_four(
    hr_variant_four: VariantCoordinates,
    hr_five_utr_sequence: str,
    hr_five_utr: FiveUTRCoordinates,        
):
    vc = prepare_alt_seq(hr_variant_four, hr_five_utr_sequence, hr_five_utr)

    assert hr_five_utr_sequence[302] == "A"
    assert vc[302] == "G"