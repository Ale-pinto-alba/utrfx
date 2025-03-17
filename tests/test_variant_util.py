import os
import pytest

from utrfx.variant_util import PrepareAltSeq, VCFfile
from utrfx.model import FiveUTRCoordinates
from utrfx.genome import Contig, GenomicRegion, Strand, VariantCoordinates, GenomeBuild

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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
            assert actual == expected

        def test_check_variant_in_cdna(
            self,
            cdna_seq_positive: str,
            five_utr_coordinates_positive: FiveUTRCoordinates,
        ):
            vc = TestPrepareAltSeq.make_variant(100, "C", "T")
            vc_instance = PrepareAltSeq(vc, cdna_seq_positive, five_utr_coordinates_positive)
            vc_pos = vc_instance.variant_position()
            assert vc_instance.check_variant_in_cdna(vc_pos) == "Variant not in the 5'UTR of the given transcript"

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
            vc_instance = PrepareAltSeq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
            vc_instance = PrepareAltSeq(vc, cdna_seq_negative, five_utr_coordinates_negative)
            vc_pos = vc_instance.variant_position()
            actual = vc_instance.prepare_alt_seq(vc_pos)
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
    

class TestVCFFile:

    @pytest.fixture(scope="class")
    def vcf_fpath(self, fpath_data_dir: str) -> str:
        return os.path.join(fpath_data_dir, "gnomad.genomes.v4.1.sites.chr8.sample.vcf.gz")

    @pytest.fixture(scope="class")
    def contig(self, genome_build: GenomeBuild) -> Contig:
        return genome_build.contig_by_name("8")

    def test_retrieve_variants_of_region(
        self,
        contig: Contig,
        vcf_fpath: str,
    ):
        # contig = genome_build.contig_by_name("8")
        with VCFfile(vcf_fpath) as vcf_fh:
            variants = vcf_fh.retrieve_variants_of_region(contig=contig, start=22_130_651, end=22_130_692)

        assert len(variants) == 9
        
        # We're getting a collection but we'd like to check the first and last variant.
        # Therefore, let's wrap the collection into a tuple to simplify testing.
        variants = tuple(variants)

        first = variants[0]
        assert first.start == 22_130_651
        assert first.end == 22_130_652

        last = variants[-1]
        assert last.start == 22_130_691
        assert last.end == 22_130_692

    def test_get_allele_frequency(
        self,
        contig: Contig,
        vcf_fpath: str,
    ):
        with VCFfile(vcf_fpath) as vcf_fh:

            af = vcf_fh.get_allele_frequency(VariantCoordinates.from_vcf_literal(contig=contig, pos=22_130_611, ref="C", alt="T"))     
            assert af == pytest.approx(1.31e-05, abs=1e-7)

            af = vcf_fh.get_allele_frequency(VariantCoordinates.from_vcf_literal(contig=contig, pos=22_130_655, ref="C", alt="T"))
            assert af == None

    def test_raises_if_not_used_as_a_context_manager(
        self,
        genome_build: GenomeBuild,
        vcf_fpath: str,
    ):
        vcf = VCFfile(vcf_fpath=vcf_fpath)
        
        contig = genome_build.contig_by_name("8")
        assert contig is not None

        with pytest.raises(AssertionError) as e:
            _ = vcf.retrieve_variants_of_region(contig=contig, start=22_130_651, end=22_130_692)
            
        assert e.value.args == ("VCFfile must be used as a context manager",)
        