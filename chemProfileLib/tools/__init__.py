from . import cutadapt, starbin, samtools, intersect

tool_list = {
    "cutadapt": cutadapt,
    "star": starbin,
    "samtools": samtools,
    "bedtools": intersect,
}