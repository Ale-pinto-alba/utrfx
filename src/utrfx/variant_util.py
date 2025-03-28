import typing
import os

import pysam

from utrfx.genome import VariantCoordinates, Strand, Contig
from utrfx.model import FiveUTRCoordinates

class AltAlleleSeq:
    """
    `AltAlleleSeq` allows the obtention of the cDNA sequence of a variant.

    It permits to check if the variant is located in the given transcript's 5'UTR.
    """
    def __init__(
        self,
        variant: VariantCoordinates,
        cdna: str,
        five_utrs: FiveUTRCoordinates,
    ):
        self._variant = variant
        self._cdna = cdna
        self._five_utrs = five_utrs
        self._gene_strand = self._obtain_gene_strand()
        self._ref = self._obtain_corresponding_ref_allele()
        self._alt = self._obtain_corresponding_alt_allele()
        self._variant_region_start = self._obtain_corresponding_variant_start()
        self._variant_cdna_pos = self.variant_position()

    def variant_position(self) -> int:
        """
        Obtain the variant position index whithin the cDNA transcript sequence as a integer.
        """
        five_utrs_tuple = []
        for region in self._five_utrs.regions:
            five_utrs_tuple.append((region.start, region.end))

        five_utrs_tuples_sorted = sorted(five_utrs_tuple)
        cdna_pos = 0
        variant_cdna_pos = None
        if self._gene_strand == Strand.POSITIVE:
            for start, end in five_utrs_tuples_sorted:
                five_utr_region_length = end - start

                if start <= self._variant_region_start <= end:
                    variant_cdna_pos = cdna_pos + (self._variant_region_start - start)
                cdna_pos += five_utr_region_length
            
            if variant_cdna_pos is not None:
                return variant_cdna_pos
        else:
            five_utrs_tuples_sorted_reversed = five_utrs_tuples_sorted[::-1]
            for start, end in five_utrs_tuples_sorted_reversed:
                five_utr_region_length = end - start
                
                if start <= self._variant_region_start <= end:
                    variant_cdna_pos = cdna_pos + (end - self._variant_region_start)
                cdna_pos += five_utr_region_length
            
            if variant_cdna_pos is not None:
                relative_variant_position = len(self._cdna) - variant_cdna_pos
                return relative_variant_position
                
    def prepare_alt_seq(self) -> str:
        """
        Get the 5'UTR region with an alternative allele from the reference cDNA.
        """
        return self._cdna[:self._variant_cdna_pos] + self._alt + self._cdna[self._variant_cdna_pos + len(self._ref):]
    
    def _obtain_gene_strand(self) -> Strand:
        for region in self._five_utrs.regions:
            return region.strand

    def _obtain_corresponding_ref_allele(self) -> str:
        if self._variant.region.strand == self._gene_strand:
            return self._variant.ref
        else:
            return reverse_complement(self._variant.ref)
        
    def _obtain_corresponding_alt_allele(self) -> str:
        if self._variant.region.strand == self._gene_strand:
            return self._variant.alt
        else:
            return reverse_complement(self._variant.alt)
        
    def _obtain_corresponding_variant_start(self) -> str:
        if self._variant.region.strand == self._gene_strand:
            return self._variant.start
        else:
            return self._variant.region.start_on_strand(self._gene_strand)
        
    def check_variant_in_cdna(self) -> str:
        """
        Check if the variant is in the 5'UTR and if the reference alleles match.
        """
        for region in self._five_utrs.regions:
            five_utr_contig = region.contig
            gene_strand = region.strand

        if self._variant.region.strand == gene_strand:
            ref = self._variant.ref
        else:
            ref = reverse_complement(self._variant.ref)

        in_variant = any(region.overlaps_with(self._variant.region) for region in self._five_utrs.regions)
        if five_utr_contig != self._variant.region.contig:
            return "Variant and cDNA not in the same contig"
        else:
            if in_variant is not True:
                return "Variant not in the 5'UTR of the given transcript"
            else:
                if self._cdna[self._variant_cdna_pos:self._variant_cdna_pos + len(ref)] != ref:
                    return "Reference alleles do not match"
                
def reverse_complement(seq: str,) -> str:
    return seq.translate(str.maketrans("ATCG", "TAGC"))


class VCFfile:
    """
    `VCFfile` represents a VCF file and allow to search for specific variants within it.

    The object must be used as a context manager to ensure proper resource cleanup.
    """
    def __init__(
        self,
        vcf_fpath: str,
    ):
        self._vcf_fpath, self._vcf_tbi_fpath = VCFfile._check_path(vcf_fpath)
        self._vcf_file = None

    @staticmethod
    def _check_path(fpath) -> typing.Tuple[str, str]:
        if os.path.isfile(fpath):
            index_fpath = fpath + ".tbi"
            if os.path.isfile(index_fpath):
                return fpath, index_fpath
            else:
                raise ValueError(f"`{index_fpath}` is not a file")
        else:
            raise ValueError(f"`{fpath}` is not a file")

    def _open_vcf(self):
        return pysam.VariantFile(
            self._vcf_fpath,
            mode="r",
            index_filename=self._vcf_tbi_fpath,
        )
    
    def __enter__(self) -> "VCFfile":
        self._vcf_file = self._open_vcf()
        return self
    
    def __exit__(self, _exc_type, _exc_value, _exc_traceback):
        self._vcf_file.close()
        self._vcf_file = None

    def retrieve_variants_of_region(
        self, 
        contig: Contig, 
        start: int,
        end: int,
    ) -> typing.Collection[VariantCoordinates]:
        """
        Get variant for the query region.

        :param contig: the query region contig.
        :param start: 0-based (excluded) start coordinate of the query region.
        :param start: 0-based (included) end coordinate of the query region.
        """
        assert self._vcf_file is not None, "VCFfile must be used as a context manager"
        variant_list = []
        for rec in self._vcf_file.fetch(contig.ucsc_name, start, end):
            pos = rec.pos
            ref = rec.ref
            for alt in rec.alts:
                variant_list.append(VariantCoordinates.from_vcf_literal(contig=contig, pos=pos, ref=ref, alt=alt))
        return variant_list
    
    def get_allele_frequency(
        self,
        variant: VariantCoordinates,
    ) -> typing.Optional[float]:
        """
        Get the allele frequency for the specified variant.

        :param variant: single variant as `VariantCoordinates`class instance.
        """
        assert self._vcf_file is not None, "VCFfile must be used as a context manager"
        contig = f"chr{variant.chrom}"
        for rec in self._vcf_file.fetch(contig, variant.start, variant.end):
            af_tuple = rec.info.get('AF', None)
            if len(af_tuple) == 1:
                return af_tuple[0] 
            else:
                af_index = - 1
                for alt in rec.alts:
                    af_index += 1
                    if alt == variant.alt:
                        af_tuple[af_index]
        return None