import os
import time

from chemProfileLib import Configuration
from chemProfileLib.tools import cutadapt, compress, starbin, fix5prime, intersect, modcount

def align(configuration: "Path to the configuration file.",
          force_overwrite: "Set to True if you want to overwrite existing files" = False,
          unstopped_keep: "Removes Stop Adapters" = False):
    """
    Takes a configured file and searches for fastq files within the set input directory. It then does the following
    procedures: Cutting adapters, aligning to the genome and intersecting the matches with a gene information file.
    As a final step, it takes the counted
    """
    config = Configuration(configuration)
    output_path = config.extend_output_path("")
    config.force_overwrite = force_overwrite

    total_time = time.time()

    if os.path.exists(output_path) == False:
        os.mkdir(output_path)

    for sample in config.samples:
        #if sample.split(".")[0] != "Y2":
        #    continue

        sample_stime = time.time()

        sample = config.extend_input_path(sample)
        print("Running sample {}".format(sample))

        # Check if compressed or not
        if sample.endswith(".gz") is False:
            sample = compress.run(config, sample)
        else:
            print("No compression needed.")

        # Cut adapters
        outfile = cutadapt.run(config, sample, remove_stopadapters=not unstopped_keep)
        # Align and convert directly to sorted bam
        outfile = starbin.run(config, outfile, bam_direct = True)
        # 5p fix (it looks like I don't need this? Nothing is fixed...)
        # outfile = fix5prime.run(config, outfile)
        # Intersect to add gene annotation information
        outfile = intersect.run(config, outfile)
        # Count modifications etcetera
        outfile = modcount.run(config, outfile)

        sample_etime = time.time()

        print("\tSample Running time: {}".format(sample_etime - sample_stime))
        print()

    end_time = time.time()

    print("Total running time: {}".format(end_time - total_time))