import typing
import os

import pysam

from utrfx.genome import VariantCoordinates, Strand, Contig, Region
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
                variant = VariantCoordinates.from_vcf_literal(contig=contig, pos=pos, ref=ref, alt=alt)
                variant.change_length = len(ref) - len(alt)
                variant_list.append(variant)
        return variant_list
      
    def get_allele_frequency(
        self,
        variant: VariantCoordinates,
    ) -> typing.Optional[float]:
        """
        Get the allele frequency for the specified variant.

        :param variant: single variant as `VariantCoordinates`class instance.
        """
        assert self._vcf_file is not None, "VCF file must be used as a context manager"

        af_field = self._get_af_field()
        if af_field is None:
            return None 

        contig = f"chr{variant.chrom}"
        for rec in self._vcf_file.fetch(contig, variant.start, variant.end):
            af_tuple = rec.info.get(af_field, None)
            if af_tuple is None:
                continue
            for af_index, alt in enumerate(rec.alts):
                if alt == variant.alt:
                    return af_tuple[af_index]
        return None
    
    def _get_af_field(self):
        """
        Determines which AF field is present in the VCF header: 'AF_joint' or 'AF'.
        """
        assert self._vcf_file is not None, "VCF file must be used as a context manager"
        if 'AF_joint' in self._vcf_file.header.info:
            return 'AF_joint'
        elif 'AF' in self._vcf_file.header.info:
            return 'AF'
        else:
            return None

class VariantClassifier:
    """
    `MutationClassifier` allows the determination of the mutation type of a given variant and which uORF is affected.
    """
    def __init__(
        self,
        canonical_uorfs_coordinates_list: typing.Collection[Region],
        canonical_uorfs_lengths_list: typing.Collection[int],
        variant_uorfs_lengths_list: typing.Collection[int],
        canonical_uorfs_ouorf_list: typing.Collection[bool],
        variant_uorfs_ouorf_list: typing.Collection[bool],
        uorf_end_pos_list: typing.Collection[int],
        variant_cdna_pos: int,
    ):
        """
        :param canonical_uorfs_lengths_list: list containing the uORFs lengths of the canonical sequence.
        :param variant_uorfs_lengths_list: list containing the uORFs lengths of the variant sequence.
        :param canonical_uorfs_ouorf_list: list containing if the uORFs of the canonical sequence are overlapping.
        :param variant_uorfs_ouorf_list: list containing if the uORFs of the variant sequence are overlapping.
        :param uorf_end_pos: integer corresponding to the end of the uORF.
        :param variant_cdna_pos: integer corresponding to the variant position within the cDNA sequence.
        """
        self._canonical_uorf_coordinates_list = canonical_uorfs_coordinates_list
        self._canonical_uorfs_lengths_list = canonical_uorfs_lengths_list
        self._variant_uorfs_lengths_list = variant_uorfs_lengths_list
        self._canonical_uorfs_ouorf_list = canonical_uorfs_ouorf_list
        self._variant_uorfs_ouorf_list = variant_uorfs_ouorf_list
        self._uorf_end_pos_list = uorf_end_pos_list
        self._variant_cdna_pos = variant_cdna_pos

    def _check_if_mutation_in_uorf(self) -> typing.Optional[str]:
        variant_in_uorf = False
        for uorf_region in self._canonical_uorf_coordinates_list:
            if uorf_region.start <= self._variant_cdna_pos <= uorf_region.end:
                variant_in_uorf = True
        if variant_in_uorf == False:
            return "Variant does not affect any canonical uORF"

    def _mutation_classifier_start_codon(self) -> typing.Optional[str]:
        """
        Compare the number of uORFs of a variant with the one of its canonical transcript.
        """
        if len(self._canonical_uorfs_lengths_list) > len(self._variant_uorfs_lengths_list):
            return "Start codon loss mutation"
        elif len(self._canonical_uorfs_lengths_list) < len(self._variant_uorfs_lengths_list):
            return "Start codon gain mutation"
        
    def _mutation_classifier_stop_codon_loss(self) -> typing.Optional[str]:
        """
        Compare each uORF of the canonical with its corresponding one in the variant sequence, if the uORF is overlapping in the variant and not in the canonical 
        and the distance between the variant position and the uORF end is less or equal to two, it would mean that the variant causes 
        the stop codon uORF loss. 
        """
        for idx, (ouorf1, ouorf2) in enumerate(zip(self._canonical_uorfs_ouorf_list, self._variant_uorfs_ouorf_list)):
            if ouorf2 is True and ouorf1 is False:
                end_pos = self._uorf_end_pos_list[idx]
                if 0 <= (end_pos - self._variant_cdna_pos) <= 2:
                    return "Stop codon loss mutation"

    def _mutation_classifier_stop_codon_gain(self) -> typing.Optional[str]:
        """
        Check if a uORF is not overlapping and if this uORF length is shorter than the one in the canonical sequence,
        meaning that a stop codon appeared because of the variant.
        """
        for uorf_index, ouorf in enumerate(self._variant_uorfs_ouorf_list):
            if ouorf is False and self._canonical_uorfs_lengths_list[uorf_index] > self._variant_uorfs_lengths_list[uorf_index]:
                return "Stop codon gain mutation"

    def _mutation_classifier_indel_snv(self) -> str:
        """
        Check if there is any change between the lengths of uORFs and if not.
        """
        for length1, length2 in zip(self._canonical_uorfs_lengths_list, self._variant_uorfs_lengths_list):
            if length1 > length2:
                return "Deletion"
            elif length1 < length2:
                return "Insertion"
        return "SNV or MNV"

    def perform_mutation_analysis(self) -> str:
        """
        Check which kind of variant is it.
        """
        variant_not_in_uorf = self._check_if_mutation_in_uorf()
        if variant_not_in_uorf:
            return variant_not_in_uorf
        start_codon_mutation = self._mutation_classifier_start_codon()
        if start_codon_mutation:
            return start_codon_mutation
        stop_codon_loss_mutation = self._mutation_classifier_stop_codon_loss()
        if stop_codon_loss_mutation:
            return stop_codon_loss_mutation
        stop_codon_gain_mutation = self._mutation_classifier_stop_codon_gain()
        if stop_codon_gain_mutation:
            return stop_codon_gain_mutation
        indel_snv_mutation = self._mutation_classifier_indel_snv()
        if indel_snv_mutation:
            return indel_snv_mutation
        
    def uorf_affected(self) -> typing.Optional[int]:
        """
        Determine which uORF is affected by the variant.
        """
        for uorf_index, uorf_region in enumerate(self._canonical_uorf_coordinates_list):
            if uorf_region.start <= self._variant_cdna_pos <= uorf_region.end:
                return uorf_index + 1
        return None
    
class VariantAA:

    def __init__(
        self,
        five_prime_seq: str,
        uorf_coordinates: Region,
        variant_cdna_pos: int
    ):
        self._five_prime_seq = five_prime_seq
        self._uorf_coordinates = uorf_coordinates
        self._variant_cdna_pos = variant_cdna_pos
    
    def variant_codon(self) -> typing.Optional[str]:
        """
        Return the codon that contains the variant.

        Work only for SNV.
        """
        for i in range(self._uorf_coordinates.start, self._uorf_coordinates.end, 3):
            codon = self._five_prime_seq[i:i + 3]
            if i <= self._variant_cdna_pos < i + 3:
                return codon
        return None
    
    @staticmethod
    def variant_amino_acid(codon: str) -> str:
        """
        Retrieve the amino acid coded by the codon.
        """
        aa_dict = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        } 
        if codon in aa_dict.keys():
            return aa_dict[codon]
    
    @staticmethod
    def amino_acids_difference_score(
        amino_acid_one: str,
        amino_acid_two: str,
    ) -> typing.Optional[int]:
        """
        Obtain the value corresponding to the BLOSUM62 value between two amino acids.
        """
        blosum62 = {
            'A': [4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
            'C': [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
            'D': [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
            'E': [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
            'F': [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
            'G': [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
            'H': [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
            'I': [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
            'K': [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
            'L': [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
            'M': [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
            'N': [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
            'P': [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
            'Q': [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
            'R': [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
            'S': [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
            'T': [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
            'V': [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
            'W': [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
            'Y': [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7],
        }
        amino_acids = list(blosum62.keys())
        index = amino_acids.index(amino_acid_two)
        return blosum62[amino_acid_one][index]
    
    @staticmethod
    def codon_usage(codon) -> float:
        """
        Retrieve the codon usage frequency in 1000.

        Krishnamurthy Subramanian, Bryan Payne, Felix Feyertag, David Alvarez-Ponce, 
        The Codon Statistics Database: A Database of Codon Usage Bias, Molecular Biology and Evolution, Volume 39, Issue 8, August 2022, msac157, 
        https://doi.org/10.1093/molbev/msac157
        """
        codon_frequency = {
            "AAA": 25.1901, "AAG": 31.4481, "AAT": 17.0444, "AAC": 18.3501,
            "ACA": 15.3184, "ACC": 18.4242, "ACG": 5.9236, "ACT": 13.5182,
            "AGA": 12.3146, "AGC": 19.9075, "AGG": 12.1033, "AGT": 12.6984,
            "ATA": 7.5636, "ATC": 19.3821, "ATT": 15.7146, "ATG": 21.0300,

            "CAA": 12.8031, "CAG": 34.5800, "CAT": 11.1677, "CAC": 15.0600,
            "CCA": 17.5929, "CCC": 20.4863, "CCG": 7.4462, "CCT": 18.2084,
            "CGA": 6.0584, "CGC": 10.3363, "CGG": 11.3963, "CGT": 4.4912,
            "CTA": 7.1547, "CTC": 19.0061, "CTG": 38.7259, "CTT": 13.3863,

        
            "GAA": 30.4087, "GAG": 40.0424, "GAT": 22.0502, "GAC": 24.8107,
            "GCA": 16.1205, "GCC": 27.8499, "GCG": 7.6763, "GCT": 18.3994,
            "GGA": 16.6903, "GGC": 22.2317, "GGG": 16.4639, "GGT": 10.6731,
            "GTA": 7.1527, "GTC": 13.9037, "GTG": 27.0722, "GTT": 10.9735,

            "TAA": 0.4856, "TAG": 0.3822, "TAT": 11.8554, "TAC": 14.3091,
            "TCA": 12.9278, "TCC": 17.7956, "TCG": 4.5676, "TCT": 15.7052, 
            "TGA": 0.8420, "TGC": 12.3082, "TGG": 12.2384, "TGT": 10.7211, 
            "TTA": 7.8952, "TTC": 19.1588, "TTG": 12.9418, "TTT": 16.9475,
        }
        return codon_frequency[codon]