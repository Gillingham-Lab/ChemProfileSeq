import os
import subprocess
from .exceptions import ToolNotFoundException
from . import compress

def check():
    """ Checks if cutadapt is installed"""
    try:
        subprocess.check_output(["STAR", "--version"], stderr=subprocess.STDOUT)
        return True
    except FileNotFoundError:
        raise ToolNotFoundException
    except PermissionError:
        raise ToolNotFoundException

def run(config, infile, bam_direct = False):
    print("\tAligning to genome")

    sampleName =  os.path.basename(infile).split(".")[0]
    outputPrefix = config.extend_output_path("")
    outputNameLog = config.extend_output_path("{sample}.Aligned.log".format(sample=sampleName))
    outputSJStatistics = config.extend_output_path("{sample}.AlignedSJStats.tab".format(sample=sampleName))



    if bam_direct:
        outputName = config.extend_output_path("{sample}.Aligned.bam".format(sample=sampleName))

        # "--outSAMprimaryFlag AllBestScore " \

        alignment = "STAR --runThreadN {alignmentThreads:d} " \
                    "--genomeDir {genomeDir} " \
                    "--readFilesIn {infile} " \
                    "--outFileNamePrefix {outprefix} " \
                    "--outSAMmultNmax 1 --outMultimapperOrder Random " \
                    "--outSAMtype BAM SortedByCoordinate " \
                    "--outFilterIntronMotifs RemoveNoncanonicalUnannotated " \
                    "--outSJfilterCountUniqueMin 300 100 100 100 --outSJfilterCountTotalMin 300 100 100 100  " \
                    "--outSJfilterReads Unique " \
                    "--limitOutSJcollapsed 2000000".format(
            alignmentThreads=max(1, config.get("general", "threads")),
            genomeDir=config.extend_path(config.get("general", "genomeDir")),
            infile=infile,
            outprefix=outputPrefix
        )

        if config.force_overwrite or os.path.exists(outputName) is False:
            print("\t\tRunning alignment")

            config.execute(alignment)

            # Move files around
            os.rename("{}Aligned.sortedByCoord.out.bam".format(outputPrefix), outputName)
            os.rename("{}Log.final.out".format(outputPrefix), outputNameLog)
            os.rename("{}SJ.out.tab".format(outputPrefix), outputSJStatistics)
            os.unlink("{}Log.progress.out".format(outputPrefix))
            #os.unlink("{}Log.out".format(outputPrefix))
        else:
            print("\t\tSkipping alignment (file already exists)")
    else:
        outputName = config.extend_output_path("{sample}.Aligned.sam".format(sample=sampleName))

        alignment = "STAR --runThreadN {alignmentThreads:d} " \
                    "--genomeDir {genomeDir} " \
                    "--readFilesIn {infile} " \
                    "--outFileNamePrefix {outprefix} " \
                    "--outSAMmultNmax 10 " \
                    "--outFilterIntronMotifs RemoveNoncanonicalUnannotated " \
                    "--outSJfilterCountUniqueMin 300 100 100 100 --outSJfilterCountTotalMin 300 100 100 100  " \
                    "--limitOutSJcollapsed 2000000".format(
            alignmentThreads=max(1, config.get("general", "threads")),
            genomeDir=config.extend_path(config.get("general", "genomeDir")),
            infile=infile,
            outprefix=outputPrefix
        )

        if config.force_overwrite or os.path.exists(outputName) is False:
            print("\t\tRunning alignment")

            config.execute(alignment)

            # Move files around
            os.rename("{}Aligned.out.sam".format(outputPrefix), outputName)
            os.rename("{}Log.final.out".format(outputPrefix), outputNameLog)
            os.rename("{}SJ.out.tab".format(outputPrefix), outputSJStatistics)
            os.unlink("{}Log.progress.out".format(outputPrefix))
            os.unlink("{}Log.out".format(outputPrefix))
        else:
            print("\t\tSkipping alignment (file already exists)")

    return outputName