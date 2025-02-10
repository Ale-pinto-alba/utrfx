from utrfx.genome import VariantCoordinates, Strand
from utrfx.model import FiveUTRCoordinates

def prepare_alt_seq(
    variant: VariantCoordinates,
    cdna: str,
    five_utrs: FiveUTRCoordinates,
) -> str:
    five_utrs_tuples = []
    for region in five_utrs.regions:
        five_utrs_tuples.append((region.start, region.end))
        for _ in range(1):
            gene_chrom = region.contig
            gene_strand = region.strand
            break

    assert variant.region.contig == gene_chrom, "Variant and 5'UTR regions not in the same contig."
    
    variant_in_five_utr_strand = variant.region.with_strand(gene_strand)

    in_variant = any(region.overlaps_with(variant.region) for region in five_utrs.regions)

    assert in_variant is True, "Variant not in the 5'UTR of the given transcript."

    five_utrs_tuples_sorted = sorted(five_utrs_tuples, key= lambda tuple:(tuple[0], tuple[1]))

    if gene_strand == Strand.POSITIVE:

        cdna_pos = 0
    
        for start, end in five_utrs_tuples_sorted:
            
            five_utr_region_length = end - start
            
            if start <= variant_in_five_utr_strand.start <= end:
                variant_cdna_pos = cdna_pos + (variant_in_five_utr_strand.start - start)
            
            cdna_pos += five_utr_region_length
        
        # assert cdna[variant_cdna_pos] == variant.ref

        return cdna[:variant_cdna_pos] + variant.alt + cdna[variant_cdna_pos + 1:]
    
    else:

        cdna_pos = 0

        five_utrs_tuples_sorted_reversed = five_utrs_tuples_sorted[::-1]
    
        for start, end in five_utrs_tuples_sorted_reversed:
            
            five_utr_region_length = end - start
            
            if start <= variant_in_five_utr_strand.start <= end:
                variant_cdna_pos = cdna_pos + (end - variant_in_five_utr_strand.start)
            
            cdna_pos += five_utr_region_length
        
        relative_variant_position = len(cdna) - variant_cdna_pos

        # assert cdna[relative_variant_position] == variant.ref.translate(str.maketrans("ATCG", "TAGC"))

        return cdna[:relative_variant_position] + variant.alt + cdna[relative_variant_position + 1:]