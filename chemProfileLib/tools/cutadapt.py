import os, subprocess
from .exceptions import ToolNotFoundException
from . import compress

def check():
    """ Checks if cutadapt is installed"""
    try:
        subprocess.check_output(["cutadapt", "--version"], stderr=subprocess.STDOUT)
        return True
    except FileNotFoundError:
        raise ToolNotFoundException

def run(config, infile, remove_stopadapters = True):
    print("\tCutting adapters")
    sampleName =  os.path.basename(infile).split(".")[0]

    outFile = infile

    if len(config.get("adapters", "3p")) > 0:
        outFile3Prime = config.extend_output_path("{sample}.3pTrimmed.fastq".format(sample=sampleName))
        outFile3PrimeLog = config.extend_output_path("{sample}.3pTrimmed.log".format(sample=sampleName))
        outFile = outFile3Prime

        cut3prime = "cutadapt -m {minLength} -a {adapter3p} {infile} > {outfile} 2> {logfile}".format(
            minLength=config.get("cutadapt", "minReadLength"),
            adapter3p=config.get("adapters", "3p"),
            infile=infile,
            outfile=outFile3Prime,
            logfile=outFile3PrimeLog,
        )

        if config.force_overwrite or (
                os.path.exists(outFile3Prime) is False and os.path.exists(outFile3Prime + ".gz") is False):
            print("\t\tCut away 3' adapters")
            config.execute(cut3prime)
        else:
            print("\t\tFound existing file, skipping 3' adapter cut.")
    else:
        print("\t\tNo 3p adapter.")

    if len(config.get("adapters", "5p")) > 0:
        outFile5PStopped = config.extend_output_path("{sample}.5pStopped.fastq".format(sample=sampleName))
        outFile5PUnStopped = config.extend_output_path("{sample}.5pUnStopped.fastq".format(sample=sampleName))
        outFile5PStoppedLog = config.extend_output_path("{sample}.5pStopped.log".format(sample=sampleName))
        outFile = outFile5PStopped

        if remove_stopadapters:
            removeReadsWith5pAdapter = "cutadapt -g {adapter5p} --untrimmed-output {outfile} {infile} > {altfile} 2> {logfile}".format(
                adapter5p=config.get("adapters", "5p"),
                infile=outFile3Prime,
                outfile=outFile5PStopped,
                altfile=outFile5PUnStopped,
                logfile=outFile5PStoppedLog
            )
        else:
            removeReadsWith5pAdapter = "cutadapt -m {minLength} -g {adapter5p} {infile} > {outfile} 2> {logfile}".format(
                minLength=config.get("cutadapt", "minReadLength"),
                adapter5p=config.get("adapters", "5p"),
                infile=outFile3Prime,
                outfile=outFile5PStopped,
                logfile=outFile5PStoppedLog
            )

        if config.force_overwrite or (
                os.path.exists(outFile5PStopped) is False and os.path.exists(outFile5PStopped + ".gz") is False):
            print("\t\tFilter out reads with an existing 5' adapter and keeping only true stopped ones.")
            config.execute(removeReadsWith5pAdapter)
            if remove_stopadapters:
                compress.run(config, outFile5PUnStopped)
            compress.run(config, outFile3Prime)
        else:
            print("\t\tFound existing file, skipping 5' adapter filter.")
    else:
        print("\t\tNo 5p adapter.")

    # We still need to unpack it in that case
    if outFile == infile and outFile.endswith(".gz"):
        outFile = config.extend_output_path("{sample}.fastq".format(sample=sampleName))
        uncompress = "gunzip -k -c {infile} > {outFile}".format(infile=infile, outFile=outFile)

        if config.force_overwrite or os.path.exists(outFile) is False:
            config.execute(uncompress)
            print("\t\tReads uncompressed to output directory.")
        else:
            print("\t\tReads already uncompressed, skipping.")

    return outFile