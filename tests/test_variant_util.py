import os
import pytest
import pysam

from utrfx.variant_util import prepare_alt_seq, VCFfile
from utrfx.model import FiveUTRCoordinates
from utrfx.genome import Contig, GenomicRegion, Strand, VariantCoordinates, GenomeBuild, Region

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
    class TestPositiveStrand:
        @pytest.fixture(scope="class")
        def five_utr_coordinates_positive(self) -> FiveUTRCoordinates:
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
        def cdna_seq_positive(self) -> str:
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
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(20, "C", "T")

            #                    *
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCTGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_del(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(20, "CGG", "C")

            #                    ***
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCGGGTTTTTAAAAACCCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_ins(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(20, "C", "CTT")

            #                    *
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCTTGGGGGTTTTTAAAAACCCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_mnv(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(20, "CGG", "CA")

            #                    ***
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCAGGGTTTTTAAAAACCCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_snp_falls_on_second_five_utr_region(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(40, "A", "T")

            #                                   *
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCGGGGGTTTTTAAAATCCCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_del_falls_on_second_five_utr_region(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(40, "ACC", "A")

            #                                   *
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCGGGGGTTTTTAAAAACCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_ins_falls_on_second_five_utr_region(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(40, "A", "ATT")

            #                                   *
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCGGGGGTTTTTAAAAATTCCCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_mnv_falls_on_second_five_utr_region(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(40, "ACC", "AT")

            #                                   *
            #     ref:  AAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTT
            expected = "AAAAACCCCCGGGGGTTTTTAAAAATCCCGGGGGTTTTT"
            actual = prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            assert actual == expected

        def test_variant_out_of_five_utr(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(100, "C", "T")
            with pytest.raises(AssertionError) as e:
                prepare_alt_seq(vc, cdna_seq_positive, five_utr_coordinates_positive)

            assert e.value.args == ("Variant not in the 5'UTR of the given transcript.",)

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
    

@pytest.fixture(scope="session")
def vcf_fpath(fpath_data_dir: str) -> str:
    return os.path.join(fpath_data_dir, "gnomad.genomes.v4.1.sites.chr8.sample.vcf.gz")

@pytest.fixture(scope="session")
def vcf_file(vcf_fpath: str) -> VCFfile:
    return VCFfile(vcf_fpath=vcf_fpath)


def test_variants_in_region(genome_build: GenomeBuild, vcf_file: VCFfile):
    contig = genome_build.contig_by_name("8")
    variants_list = vcf_file.retrieve_variants_of_region(contig=contig, start=22_130_610, end=22_130_650)

    assert len(variants_list) == 9

    vcf_file.close_vcf()