import typing
import re

import pandas as pd
import numpy as np

from utrfx.genome import GenomeBuild, GenomicRegion, Strand
from utrfx.model import FiveUTRCoordinates, TranscriptCoordinates

def read_gtf_into_txs(fpath: str, genome_build: GenomeBuild) -> typing.Collection[TranscriptCoordinates]:
    """
    Parse a GTF file and return the available transcripts.
    """
    assert fpath.endswith(".gtf"), "Not a GTF file."
    gtf_df = pd.read_csv(fpath, sep = "\t", header = None, comment = "#")

    gtf_df.columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]

    fields = [
        "transcript_id",
    ]
    
    for field in fields:
        gtf_df[field] = gtf_df["attribute"].apply(lambda label: re.findall(rf'{field} "([^"]*)"', label)[0] if rf'{field} "' in label else '')
    
    pd.set_option("future.no_silent_downcasting", True)
    gtf_df.replace('', np.nan, inplace=True)
    gtf_df.drop(["source", "score", "frame", "attribute"], axis=1, inplace=True)
    assert list(gtf_df.columns) == ["seqname", "feature", "start", "end", "strand", "transcript_id"]
    
    utr_df = gtf_df[gtf_df["feature"] == "UTR"]
    start_codon_df = gtf_df[gtf_df["feature"] == "start_codon"]
    transcripts = []
    
    for transcript_id, group in utr_df.groupby("transcript_id"):
        contig = genome_build.contig_by_name(group["seqname"].iloc[0])
        if contig is None:
            print("No contig found.")
        else:
            start_codon = start_codon_df[start_codon_df["transcript_id"] == transcript_id]

            temp_utr_5prime_list = []

            if not start_codon.empty:
                temp_start_codon = start_codon.iloc[0]
                actual_start_codon_strand = parse_strand(temp_start_codon["strand"])
                temp_start_codon = GenomicRegion(
                    contig=contig,
                    start=int(temp_start_codon["start"]) - 1,
                    end=int(temp_start_codon["end"]),
                    strand=Strand.POSITIVE,
                ).with_strand(actual_start_codon_strand)

                for _, row in group.iterrows():
                    actual_feature_strand = parse_strand(row["strand"])
                    utr_region = GenomicRegion(
                        contig=contig,
                        start=row["start"] - 1, 
                        end=row["end"],
                        strand=Strand.POSITIVE,
                    ).with_strand(actual_feature_strand)

                    if utr_region.distance_to(temp_start_codon) >= 0:
                        temp_utr_5prime_list.append(utr_region)
            
            if temp_utr_5prime_list:
                transcripts.append(
                    TranscriptCoordinates(
                        tx_id=transcript_id,
                        five_utr=FiveUTRCoordinates(
                            regions=temp_utr_5prime_list,
                        ),
                    )
                )

    return transcripts

def parse_strand(val: str) -> Strand:
    if val == "+":
        return Strand.POSITIVE
    elif val == "-":
        return Strand.NEGATIVE
    else:
        raise ValueError(f"Unknown strand value {val}")