from utrfx.genome import VariantCoordinates, GenomicRegion, GenomeBuild, Strand, GRCh38
from utrfx.model import FiveUTRCoordinates

def prepare_alt_seq(
    variant: VariantCoordinates,
    genomic_sequence: str,
    five_utrs: FiveUTRCoordinates,
) -> str:
    five_utrs_tuples = []
    for region in five_utrs.regions:
        five_utrs_tuples.append((region.start, region.end))
        for _ in range(1):
            gene_strand = region.strand
            gene_chrom = region.contig
            break

    variant_contig = variant.chrom
    variant_genomic_region = GenomicRegion(
                                contig=GenomeBuild.contig_by_name(GRCh38, name=variant_contig),
                                start=variant.start - 1,
                                end=variant.end,
                                strand=Strand.POSITIVE
                                ).with_strand(other=gene_strand)
    
    assert variant_genomic_region.contig == gene_chrom

    five_utrs_tuples.sort()
    variant_in_five_utr = False
    for start, end in five_utrs_tuples:
        if start <= variant_genomic_region.start <= end:
            variant_in_five_utr = True

    assert variant_in_five_utr is True, "Variant not in the 5'UTR of the given transcript."

    for _ in range(1):
        for start, end in five_utrs_tuples:
            variant_relative_position = variant_genomic_region.start - start
            break
    
    assert genomic_sequence[variant_relative_position + 1: variant_relative_position + len(variant.ref) + 1] == variant.ref, "Reference do not match the position."

    upstream_variant = genomic_sequence[:variant_relative_position + 1]
    downstream_variant = genomic_sequence[variant_relative_position + len(variant.ref) + 1:]

    return upstream_variant + variant.alt + downstream_variant