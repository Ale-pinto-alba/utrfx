import typing
import os

import pysam

from utrfx.genome import VariantCoordinates, Strand, Contig
from utrfx.model import FiveUTRCoordinates

def prepare_alt_seq(
    variant: VariantCoordinates,
    cdna: str,
    five_utrs: FiveUTRCoordinates,
) -> str:
    """
    Get the 5'UTR region with an alternative allele from the reference cDNA.

    :param variant: Variant as VariantCoordinates instance.
    :param cdna: 5'UTR region as cDNA (with reference allele).
    :param five_utrs: Genomic Regions of the 5'UTR. 
    """
    five_utrs_tuples = []
    for region in five_utrs.regions:
        five_utrs_tuples.append((region.start, region.end))
        gene_chrom = region.contig
        gene_strand = region.strand

    assert variant.region.contig == gene_chrom, "Variant and 5'UTR regions not in the same contig."

    in_variant = any(region.overlaps_with(variant.region) for region in five_utrs.regions)
    assert in_variant is True, "Variant not in the 5'UTR of the given transcript."

    if variant.region.strand == gene_strand:
        ref = variant.ref
        alt = variant.alt
        variant_region_start = variant.start
    else:
        ref = reverse_complement(variant.ref)
        alt = reverse_complement(variant.alt)
        variant_region_start = variant.region.start_on_strand(gene_strand)

    five_utrs_tuples_sorted = sorted(five_utrs_tuples)
    cdna_pos = 0
    if gene_strand == Strand.POSITIVE:
        for start, end in five_utrs_tuples_sorted:
            five_utr_region_length = end - start

            if start <= variant_region_start <= end:
                variant_cdna_pos = cdna_pos + (variant_region_start - start)
            cdna_pos += five_utr_region_length
        
        assert cdna[variant_cdna_pos:variant_cdna_pos + len(ref)] == ref
        return cdna[:variant_cdna_pos] + alt + cdna[variant_cdna_pos + len(ref):]
    
    else:
        five_utrs_tuples_sorted_reversed = five_utrs_tuples_sorted[::-1]
        for start, end in five_utrs_tuples_sorted_reversed:
            five_utr_region_length = end - start
            
            if start <= variant_region_start <= end:
                variant_cdna_pos = cdna_pos + (end - variant_region_start)
            cdna_pos += five_utr_region_length
        
        relative_variant_position = len(cdna) - variant_cdna_pos
        assert cdna[relative_variant_position:relative_variant_position + len(ref)] == ref
        return cdna[:relative_variant_position] + alt + cdna[relative_variant_position + len(ref):]

def reverse_complement(seq: str) -> str:
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
        contig: Contig, 
        start: int,
        end: int,
    ) -> typing.Optional[float]:
        """
        Get the allele frequency for the specified region.

        :param contig: the query region contig.
        :param start: 0-based (exclusive) start coordinate of the query region.
        :param end: 0-based (inclusive) end coordinate of the query region.
        """
        assert self._vcf_file is not None, "VCFfile must be used as a context manager"
        for rec in self._vcf_file.fetch(contig.ucsc_name, start, end):
            af = rec.info.get('AF', None)
            if isinstance(af, list) and len(af) == 1:
                return af[0] 
            elif isinstance(af, list):
                return af  
            else:
                return af
        return None