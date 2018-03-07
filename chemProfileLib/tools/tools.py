import os
import io
import subprocess
import pandas as pd
import numpy as np
from ..configuration import Configuration

def enrich_intersectedTab(config: Configuration, bam_file: str, intersect_file: str):
    # Read in the whole bam file
    tmp_file = intersect_file + "-tmp"

    """intersect_data = pd.read_csv(intersect_file, sep="\t", header=None,
        names=["chr", "read_start", "read_end", "qname", "strand", "type", "gene_start", "gene_stop", "gene_id"],
        index_col=[3, 8]
    )"""

    #intersect_data = pd.read_csv(io.BytesIO(subprocess.check_output("cat {} | head -1000".format(intersect_file), shell=True)),
    intersect_data = pd.read_csv(intersect_file,
         sep="\t", header=None,
         names=["chr", "read_start", "read_end", "qname", "strand", "type", "gene_start",
                "gene_stop", "gene_id"],
         dtype={"chr": str, "read_start": np.int32, "read_end": np.int32, "qname": str, "strand": str,
             "type": str, "gene_start": np.int32, "gene_stop": np.int32, "gene_id": str},
    )

    #cli = "samtools view {bamfile} | head -1000".format(bamfile=bam_file)
    cli = "samtools view {bamfile}".format(bamfile=bam_file)
    """bam_data = pd.read_csv(io.BytesIO(subprocess.check_output(cli, shell=True)), delim_whitespace=True, header=None,
        #names=["qname", "qflags", "chr", "rstart", "A", "CIGAR", "B", "C", "D", "read", "qual", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        names=["qname", "chr", "read_start", "CIGAR", "read"],
        dtype={"qname": str, "chr": str, "read_start": np.int32, "CIGAR": str, "read": str},
        usecols=[0, 2, 3, 5, 9]
    )

    def cleanup(x):
        cigar = x["CIGAR"]
        length = 0
        buf = ""
        for c in cigar:
            if c in "0123456789":
                buf += c
                continue

            if c in "MDNX=":
                length += int(buf)

        x["read_start"] -= 1
        x["read_end"] = x["read_start"] + length
        return x


    print("\t\t\tPrepare bam/sam file by correcting read_start to 0-indexed and adding read_end.")
    bam_data["read_end"] = 0
    bam_data = bam_data.apply(cleanup, axis=1)
    bam_data.set_index(["qname", "chr", "read_start", "read_end"], inplace=True)"""

    bam_data = pd.read_csv(io.BytesIO(subprocess.check_output(cli, shell=True)), delim_whitespace=True, header=None,
        # names=["qname", "qflags", "chr", "rstart", "A", "CIGAR", "B", "C", "D", "read", "qual", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        names=["qname", "CIGAR", "read"],
        dtype={"qname": str, "CIGAR": str, "read": str},
        usecols=[0, 5, 9],
        index_col="qname"
    )

    print(bam_data.head(2))

    #bam_data["strand"] = ""
    #bam_data.loc[bam_data["qflags"] & 0x10, "strand"] = "-"
    #bam_data.loc[bam_data["qflags"] & 0x10, "strand"] = "+"

    #print(bam_data.head(2))

    #bam_data.set_index(["qname", "chr", "read_start"], inplace=True)

    print("\t\t\tJoin Intersection with bam file data")
    #intersected = intersect_data.join(bam_data, on=["qname", "chr", "read_start", "read_end"])
    intersected = intersect_data.join(bam_data, on=["qname"])

    # Clean up old data
    del intersect_data, bam_data

    print("\t\tReading in gene ids with gene name association and joining with intersection data")
    id2name = pd.read_csv(config.extend_path(config.get("general", "idnamemap")),
        sep="\t",
        header=None, names=["gene_id", "gene_name"], index_col=["gene_id"])
    intersected = intersected.join(id2name, on="gene_id")

    # Change column order
    intersected = intersected[["gene_id", "chr", "read_start", "read_end", "qname", "strand", "type", "gene_start", "gene_stop", "CIGAR", "read", "gene_name"]]

    # Drop duplicates
    #print("\t\tTry deduplication.")
    #print("\t\t\tBefore deduplication: {}".format(len(intersected)))
    #intersected.drop_duplicates(inplace=True)
    #print("\t\t\tAfter deduplication: {}".format(len(intersected)))

    # Sort
    print("\t\tSorting the table for chr, read_start and read_end")
    intersected.sort_values(["chr", "gene_id", "read_start", "read_end"], inplace=True)

    print("\t\tWriting file.")
    intersected.to_csv(tmp_file, header=None, sep="\t", index=False)

    os.unlink(intersect_file)
    os.rename(tmp_file, intersect_file)


def countUpEverything(config: Configuration, bam_file: str):
    with open(bam_file, "r") as fo:
        for line in bam_file:
            line_parts = line.strip().split("\t")