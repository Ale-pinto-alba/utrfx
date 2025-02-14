from utrfx.genome import VariantCoordinates, Strand
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