import os
import pickle
import pysam


from ..configuration import Configuration
from ..helpers.intersection_file_processor import pool_processor


def get_output_file(config, sample):
    sampleName = os.path.basename(sample).split(".")[0]
    outFile = config.extend_output_path("{}.Modcount.data".format(sampleName))
    return outFile


def run(config: Configuration, infile):
    sampleName = os.path.basename(infile).split(".")[0]
    outFile = config.extend_output_path("{}.Modcount.data".format(sampleName))

    print("\tStarting to count stops and mutations...")

    if os.path.exists(outFile) is True and config.force_overwrite is False:
        print("\t\tModcount file found; skipping.")
        return outFile

    # Set path to genome file (index must be present)
    genome_file = config.extend_path(
        os.path.join(
            config.get("general", "genomeDir"),
            config.get("general", "genomeFile")
        )
    )

    # prepare genome
    genome = pysam.FastaFile(genome_file)

    pool = pool_processor(infile, genome, config.get("modcount", "minReadsPerGene"), int(config.get("general", "threads")) - 1)

    genes = 0
    genes_skipped = 0
    with open(outFile, "wb") as fh:
        for result in pool:
            if result is None:
                continue

            genes+=1
            # Filter a second time
            if result[1][:,1].mean() < 0.5:
                genes_skipped += 1
                continue

            pickle.dump(result, fh)

    print("\t\tFinished processing genes. Out of {} genes, {} have been saved.".format(genes, genes-genes_skipped))

    return outFile