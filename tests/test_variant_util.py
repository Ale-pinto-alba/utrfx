import pytest

# from utrfx.util import fetch_genomic_sequence_from_ensembl
from utrfx.variant_util import prepare_alt_seq
from utrfx.model import FiveUTRCoordinates
from utrfx.genome._variant import VariantCoordinates

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