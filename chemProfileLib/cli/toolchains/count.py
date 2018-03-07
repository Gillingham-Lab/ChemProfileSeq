import os
import time

from chemProfileLib import Configuration
from chemProfileLib.tools import intersect, modcount

def count(configuration: "Path to the configuration file.",
          force_overwrite: "Set to True if you want to overwrite existing files" = False):
    """
    .
    """
    config = Configuration(configuration)
    output_path = config.extend_output_path("")
    config.force_overwrite = force_overwrite

    total_time = time.time()

    if os.path.exists(output_path) == False:
        os.mkdir(output_path)

    for sample in config.samples:
        # Count modifications etcetera
        sample_stime = time.time()

        infile = config.extend_input_path(sample)
        print("Running sample {}".format(sample))

        if os.path.exists(infile):
            infile = intersect.get_output_file(config, infile)
            outfile = modcount.run(config, infile)
        else:
            print("Could not run sample (Intersection file not found, {})".format(infile))

        sample_etime = time.time()

        print("\tSample Running time: {}".format(sample_etime - sample_stime))
        print()

    end_time = time.time()

    print("Total running time: {}".format(end_time - total_time))