import enum

from ._genome import Contig, GenomicRegion, Strand

class VariantClass(enum.Enum):
    """
    `VariantClass` represents a high-level variant category
    which mostly corresponds to the structural variant categories
    of the Variant Call Format specification,
    but includes type for single nucleotide variants (SNV) and multi-nucleotide variant (MNV).
    """

    DEL = 0
    """
    A deletion - a variant with a net loss of sequence from the alternative allele
    regardless of its size.
    
    Both a deletion of 1 base pair and a deletion of 1000 base pairs are acceptable.
    """

    DUP = 1
    """
    Duplication (tandem or interspersed).
    """

    INS = 2
    """
    Insertion of a novel sequence.
    """

    INV = 3
    """
    Inversion of a chromosome segment.
    """

    MNV = 4
    """
    Multi-nucleotide variant (e.g. `AG>CT`) that is not a duplication, deletion, or insertion.
    May be called INDEL.
    """

    SNV = 5
    """
    Single nucleotide variant.
    """

    TRANSLOCATION = 6
    """
    A chromosomal translocation, which occurs when a chromosome breaks
    and the (typically two) fragmented pieces re-attach to different chromosomes.
    """


class VariantCoordinates:
    """
    A representation of coordinates of sequence and symbolic variants.

    Note, the breakend variants are not currently supported.
    """

    @staticmethod
    def from_vcf_literal(
        contig: Contig,
        pos: int,
        ref: str,
        alt: str,
    ):
        """
        Create `VariantCoordinates` from a variant in VCF literal notation.

        Note, this function must *not* be used to create a VCF symbolic variant
        (e.g. `<DEL>` or translocation).
        Use :func:`from_vcf_symbolic` instead.

        **Example**

        Create a variant from a VCF line:
        ```
        #CHROM  POS     ID  REF ALT ...
        chr1    1001    .   C   T
        ```

        We first must decide on the genome build. Most of the time, we should use GRCh38:

        >>> from utrfx.genome import GRCh38
        >>> build = GRCh38

        Then, we access the contig for ``'chr1'``:
;
        >>> chr1 = build.contig_by_name('chr1')

        Last, we create the variant coordinates:

        >>> from utrfx.genome import VariantCoordinates
        >>> vc = VariantCoordinates.from_vcf_literal(
        ...     contig=chr1, pos=1001, ref='C', alt='T',
        ... )

        Now can test the properties:

        >>> vc.start, vc.end, vc.ref, vc.alt, len(vc), vc.change_length
        (1000, 1001, 'C', 'T', 1, 0)

        Args:
            contig: a :class:`~utrfx.genome.Contig` for the chromosome
            pos: a 1-based coordinate of the first base of the reference allele, as described in VCF standard
            ref: a `str` with the REF allele. Should meet the requirements of the VCF standard.
            alt: a `str` with the ALT allele. Should meet the requirements of the VCF standard.
        """
        region = GenomicRegion(
            contig=contig,
            start=pos - 1,
            end=pos + len(ref) - 1,
            strand=Strand.POSITIVE,
        )

        change_length = len(ref) - len(alt)

        return VariantCoordinates(
            region=region,
            ref=ref,
            alt=alt,
            change_length=change_length,
        )

    @staticmethod
    def from_vcf_symbolic(
        contig: Contig,
        pos: int,
        end: int,
        ref: str,
        alt: str,
        svlen: int,
    ):
        """
        Create `VariantCoordinates` from a variant in VCF symbolic notation.

        Note, this function must *not* be used to create a VCF sequence/literal variant.
        Use :func:`from_vcf_literal` instead.

        **Example**

        Let's create a symbolic variant from a VCF line:

        ```
        #CHROM   POS      ID   REF   ALT     QUAL   FILTER   INFO
        2        321682   .    T     <DEL>   6      PASS     SVTYPE=DEL;END=321887;SVLEN=-205
        ```

        We first must decide on the genome build. Most of the time, we should use GRCh38:

        >>> from utrfx.genome import GRCh38
        >>> contig = GRCh38.contig_by_name('2')

        Now, we create the coordinates as:

        >>> vc = VariantCoordinates.from_vcf_symbolic(
        ...     contig=contig, pos=321682, end=321887,
        ...     ref='T', alt='<DEL>', svlen=-205,
        ... )

        Now can test the properties:

        >>> vc.start, vc.end, vc.ref, vc.alt, len(vc), vc.change_length
        (321681, 321887, 'T', '<DEL>', 206, -205)

        Args:
            contig: a :class:`~utrfx.genome.Contig` for the chromosome
            pos: a 1-based coordinate of the first base of the affected reference allele region
            end: a 1-based coordinate of the last base of the affected reference allele region
            ref: a `str` with the REF allele. Most of the time, it is one of `{'N', 'A', 'C', 'G', 'T'}`
            alt: a `str` with the ALT allele, e.g. one of ``{'<DEL>', '<DUP>', '<INS>', '<INV>'}``
            svlen: an `int` with change length (the difference between ref and alt allele lengths)
        """
        assert alt.startswith("<") and alt.endswith(">")
        region = GenomicRegion(
            contig=contig,
            start=pos - 1,  # convert to 0-based coordinates,
            end=end,
            strand=Strand.POSITIVE,
        )

        return VariantCoordinates(
            region=region,
            ref=ref,
            alt=alt,
            change_length=svlen,
        )

    def __init__(
        self, region: GenomicRegion,
        ref: str,
        alt: str,
        change_length: int,
    ):
        assert isinstance(region, GenomicRegion)
        self._region = region
        assert isinstance(ref, str)
        self._ref = ref
        assert isinstance(alt, str)
        self._alt = alt
        assert isinstance(change_length, int)
        self._change_length = change_length

    @property
    def chrom(self) -> str:
        """
        Get the label of the chromosome/contig where the variant is located.
        """
        return self._region.contig.name

    @property
    def start(self) -> int:
        """
        Get the 0-based start coordinate (excluded) of the first base of the :attr:`ref` allele.
        """
        return self._region.start

    @property
    def end(self) -> int:
        """
        Get the 0-based end coordinate (included) of the last base of the :attr:`ref` allele.
        """
        return self._region.end

    @property
    def region(self) -> GenomicRegion:
        """
        Get the genomic region spanned by the :attr:`ref` allele.
        """
        return self._region

    @property
    def ref(self) -> str:
        """
        Get the reference allele (e.g. "A", "CCT", "N"). The allele may be an empty string.
        """
        return self._ref

    @property
    def alt(self) -> str:
        """
        Get the alternate allele (e.g. "A", "GG", "<DEL>").

        The allele may be an empty string for sequence variants.
        The symbolic alternate allele follow the VCF notation and use the `<` and `>` characters
        (e.g. "<DEL>", "<INS:ME:SINE>").
        """
        return self._alt

    @property
    def change_length(self) -> int:
        """
        Get the change of length between the `ref` and `alt` alleles due to the variant presence.
        """
        return self._change_length

    @property
    def variant_class(self) -> VariantClass:
        """
        Get a :class:`VariantClass` category.
        """
        if self.is_structural():
            # Expecting a `str` like <DEL>, <INS>, <DUP>, <INV>, ...
            return VariantClass[self.alt[1:-1]]
        else:
            if len(self.ref) > len(self.alt):
                if self.alt == self.ref[:len(self.alt)]:
                    # alt is prefix of ref, hence a DEL
                    return VariantClass.DEL
                else:
                    return VariantClass.MNV
            elif len(self.ref) < len(self.alt):
                
                if self.ref == self.alt[:len(self.ref)]:
                    # ref is prefix of alt, hence a INS.
                    # However, it may as well be a duplication,
                    # but it's hard to say from the information on hand.
                    return VariantClass.INS
                else:
                    return VariantClass.MNV
            else:
                if len(self.ref) == 1:
                    return VariantClass.SNV
                else:
                    return VariantClass.MNV

    def is_structural(self) -> bool:
        """
        Checks if the variant coordinates use structural variant notation as described by Variant Call Format
        (`VCF <https://en.wikipedia.org/wiki/Variant_Call_Format>`_).

        Ane example of *structural* variant notation::

          chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10


        as opposed to the *sequence* (literal) notation::

          chr5  101 . NACGTACGTAC N

        :return: `True` if the variant coordinates use structural variant notation.
        """
        return (
            len(self._alt) != 0
            and self._alt.startswith("<")
            and self._alt.endswith(">")
        )

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return len(self._region)

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, VariantCoordinates)
            and self.region == other.region
            and self.ref == other.ref
            and self.alt == other.alt
            and self.change_length == other.change_length
        )

    def __hash__(self) -> int:
        return hash((self._region, self._ref, self._alt, self._change_length))

    def __str__(self) -> str:
        return (
            f"VariantCoordinates(region={self.region}, "
            f"ref={self.ref}, alt={self.alt}, "
            f"change_length={self.change_length})"
        )

    def __repr__(self) -> str:
        return str(self)
