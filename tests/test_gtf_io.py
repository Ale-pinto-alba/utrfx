import os

import pytest

from utrfx.genome import GenomeBuild, Strand
from utrfx.gtf_io import GTFio

class TestGtfIo:

    class TestFiveUTRnotexplicit:

        @pytest.fixture(scope="class")
        def fpath_example_gtf(self, fpath_data_dir: str) -> str:  
            return os.path.join(fpath_data_dir,  "Homo.sapiens.GRCh38_sample_chr22.gtf")

        def test_read_gtf_into_txs(
            self,
            fpath_example_gtf: str,
            genome_build: GenomeBuild,
        ):
            gtf_file = GTFio(fpath= fpath_example_gtf)
            transcripts = gtf_file.extract_five_utrs_if_not_explicit(genome_build=genome_build)

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
            assert len(our_another_favorite_tx.five_utr.regions) == 2
            three, four = sorted(our_another_favorite_tx.five_utr.regions, key=lambda region: region.start)

            assert three.contig.name == "22"
            assert three.start_on_strand(Strand.POSITIVE) == 23_980_411
            assert three.end_on_strand(Strand.POSITIVE) == 23_980_469
            assert three.strand == Strand.NEGATIVE
            
            assert four.contig.name == "22"
            assert four.start_on_strand(Strand.POSITIVE) == 23_976_892
            assert four.end_on_strand(Strand.POSITIVE) == 23_977_067
            assert four.strand == Strand.NEGATIVE

    class TestFiveUTRexplicit:
        
        @pytest.fixture(scope="class")
        def fpath_example_gtf(self, fpath_data_dir: str) -> str:  
            return os.path.join(fpath_data_dir,  "Homo_sapiens.GRCh38.113_sample_chr8.gtf")
        
        def test_read_gtf_into_txs(
            self,
            fpath_example_gtf: str,
            genome_build: GenomeBuild,
        ):
            gtf_file = GTFio(fpath= fpath_example_gtf)
            transcripts = gtf_file.extract_five_utrs_if_explicit(genome_build=genome_build)

            assert len(transcripts) == 3_591

            # Positive strand
            for tx in transcripts:
                if tx.tx_id == "ENST00000517969":
                    our_favorite_tx = tx
                    break

            assert our_favorite_tx is not None
            assert len(our_favorite_tx.five_utr.regions) == 2
            one, two = sorted(our_favorite_tx.five_utr.regions, key=lambda region: region.start)

            assert one.contig.name == "8"
            assert one.start == 42_154_130
            assert one.end == 42_154_536
            assert one.strand == Strand.POSITIVE

            assert two.contig.name == "8"
            assert two.start == 42_154_615
            assert two.end == 42_154_687
            assert two.strand == Strand.POSITIVE


            # Negative strand
            for tx in transcripts:
                if tx.tx_id == "ENST00000381418":
                    our_another_favorite_tx = tx
                    break

            assert our_another_favorite_tx is not None
            assert len(our_another_favorite_tx.five_utr.regions) == 2
            three, four = sorted(our_another_favorite_tx.five_utr.regions, key=lambda region: region.start)

            assert three.contig.name == "8"
            assert three.start_on_strand(Strand.POSITIVE) == 22_130_427
            assert three.end_on_strand(Strand.POSITIVE) == 22_131_010
            assert three.strand == Strand.NEGATIVE
            
            assert four.contig.name == "8"
            assert four.start_on_strand(Strand.POSITIVE) == 22_129_170
            assert four.end_on_strand(Strand.POSITIVE) == 22_129_210
            assert four.strand == Strand.NEGATIVE