import os
import subprocess

from .exceptions import ToolNotFoundException
from .tools import enrich_intersectedTab
from ..configuration import Configuration

def check():
    """ Checks if bedtools is installed """
    try:
        subprocess.check_output(["bedtools", "--version"], stderr=subprocess.STDOUT)
        return True
    except FileNotFoundError:
        raise ToolNotFoundException

def get_output_file(config, sample):
    sampleName = os.path.basename(sample).split(".")[0]
    outFile = config.extend_output_path("{}.Intersected.tab".format(sampleName))
    return outFile

def run(config: Configuration, infile: str):
    print("\tIntersecting reads with annotation data to derive gene information")
    sampleName = os.path.basename(infile).split(".")[0]

    annotationFile = config.extend_path(config.get("general", "annotationFile"))
    annotationFileType = annotationFile.split(".")[-1]
    outfile = config.extend_output_path("{sample}.Intersected.tab".format(sample=sampleName))

    awkParser = {
        "gff": """awk 'BEGIN { OFS = "\t";} {split($21,a,";"); sub(/ID=/,"",a[1]); print $1, $2, $3, $4, $6, $15, $16, $17, a[1]}'""",
        "gff3": """awk 'BEGIN { OFS = "\t";} {split($21,a,";"); sub(/ID=/,"",a[1]); print $1, $2, $3, $4, $6, $15, $16, $17, a[1]}'""",
        "gft": """awk 'BEGIN { OFS = "\t";} {id = ""; name = ""; for(i = 21; i < 40; ++i) {if ($i == "gene_id") {id = substr($(i+1), 2, length($(i+1))-3)}; if ($i == "gene_name") {name = substr($(i+1), 2, length($(i+1))-3)};}; print $1, $2, $3, $4, $6, $15, $16, $17, id "_" name;}'""",
        "bed": """awk 'BEGIN { OFS = "\t";} {print $1, $2, $3, $4, $6, "NA", $14, $15, $16}'""",
    }

    cli = "intersectBed -s -wo -split -bed -abam {infile} -b {annotations} | {awk} > {outfile}".format(
        infile=infile,
        annotations=annotationFile,
        awk=awkParser[annotationFileType],
        outfile=outfile
    )

    if config.force_overwrite or os.path.exists(outfile) is False:
        print("\t\tIntersecting...")
        config.execute(cli)
        print("\t\tAdding CIGAR and pure read")
        enrich_intersectedTab(config, infile, outfile)
    else:
        print("\t\tFile already exists, skipped.")


    return outfile
