import os
import subprocess

from .exceptions import ToolNotFoundException

def check():
    """ Checks if samtools is installed."""
    p = subprocess.Popen("samtools", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()

    if err.startswith(b"\nProgram: samtools (Tools for alignments in the SAM format)\n"):
        return True
    else:
        raise ToolNotFoundException


def samToBam(config, infile):
    print("Converting sam to bam")
    sampleName =  os.path.basename(infile).split(".")[0]
    outfile = infile[:-3] + "bam"

    cli = "samtools view -@ {threads:d} -bS {infile} > {outfile} 2> {logfile}".format(
        threads=max(1, config.get("general", "threads")),
        infile=infile,
        outfile=outfile,
        logfile="/dev/null"
    )

    config.execute(cli)

    return outfile