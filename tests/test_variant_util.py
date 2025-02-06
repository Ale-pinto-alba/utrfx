import pytest

# from utrfx.util import fetch_genomic_sequence_from_ensembl
from utrfx.variant_util import prepare_alt_seq
from utrfx.model import FiveUTRCoordinates
from utrfx.genome import Contig, GenomicRegion, Strand, VariantCoordinates


# @pytest.mark.online
# def test_adsl_variant(
#     adsl_variant: VariantCoordinates,   # (SNV) gene on + strand 
#     adsl_five_utr: FiveUTRCoordinates,
# ):
#     genomic_sequence = fetch_genomic_sequence_from_ensembl("ENST00000623063")

#     assert genomic_sequence[10] == "T"

#     variant_genomic_sequence = prepare_alt_seq(adsl_variant, genomic_sequence, adsl_five_utr)

#     assert variant_genomic_sequence[10] == "C"

#     assert len(genomic_sequence) == len(variant_genomic_sequence)
#     assert genomic_sequence[:10] == variant_genomic_sequence[:10]
#     assert genomic_sequence[11:] == variant_genomic_sequence[11:]


# @pytest.mark.online
# def test_hr_variant(
#     hr_variant: VariantCoordinates, # (SNV) gene on + strand
#     hr_five_utr: FiveUTRCoordinates,
# ):
#     genomic_sequence = fetch_genomic_sequence_from_ensembl("ENST00000381418")

#     assert genomic_sequence[405] == "A"

#     variant_genomic_sequence = prepare_alt_seq(hr_variant, genomic_sequence, hr_five_utr)

#     assert variant_genomic_sequence[405] == "G"

#     assert len(genomic_sequence) == len(variant_genomic_sequence)
#     assert genomic_sequence[:405] == variant_genomic_sequence[:405]
#     assert genomic_sequence[406:] == variant_genomic_sequence[406:]


# @pytest.mark.online
# def test_msh2_variant(
#     msh2_variant: VariantCoordinates, # Insertion
#     msh2_five_utr: FiveUTRCoordinates,
# ):
#     genomic_sequence = fetch_genomic_sequence_from_ensembl("ENST00000543555")

#     assert genomic_sequence[43] == "G"

#     variant_genomic_sequence = prepare_alt_seq(msh2_variant, genomic_sequence, msh2_five_utr)

#     assert variant_genomic_sequence[43:45] == "GA"

#     assert len(genomic_sequence) + 1 == len(variant_genomic_sequence)
#     assert genomic_sequence[:43] == variant_genomic_sequence[:43]
#     assert genomic_sequence[-10:] == variant_genomic_sequence[-10:]


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
                    start=15,
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
        """
        # Genomic coordinates (1-based):
        # 
        #      11       20        30        40        50
        #       |        |         |         |         |
        #       |                                      |
        #       |     5'UTR (1)    5'UTR (2)           |
        #       v    vvvvvvvvvv     vvvvv              v
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
