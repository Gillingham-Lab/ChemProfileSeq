import os
import csv
from typing import Tuple

from . import compress

def run(config, infile):
    sampleName =  os.path.basename(infile).split(".")[0]
    outFile = config.extend_output_path("{}.5pFixed.sam".format(sampleName))
    misMatchFile = config.extend_output_path("{}.5pMisMatches.sam".format(sampleName))

    fix_5p(infile, outFile, misMatchFile)

    compress.run(misMatchFile)

    return outFile

def fix_5p(infile, outfile, misMatchFile):
    """
    Fixes 5p mismatches by writing correcting them and writing them to outfile. If the line has been modified, the
    original line is written to misMatchFile
    :param infile: The .sam file to correct 5p mismatches
    :param outfile: Contains after the run only correct 5p matches
    :param misMatchFile: Contains the original line of lines that were fixed
    :return:
    """
    with open(infile, "r") as read, open(outfile, "w") as fixed, open(misMatchFile, "w") as mismatched:
        i = 0
        mis = 0

        for line in infile:
            originalLine = line
            lineStr, misCount = fix_line(line)

            # Only write is line is not empty
            if line:
                fixed.write(lineStr)

            # Write original line to backup file in case it was modified
            if misCount:
                mismatched.write(originalLine)
                mis += 1

            i += 1

        print("\t\tFixed {} mismatches reads (total number of reads is {})".format(mis, i))

def fix_line(line) -> Tuple[str, int]:
    # sam file format example:
    # HISEQ:108:H7N5WADXX:1:1107:18207:18800  0       gi|207113128|ref|NR_002819.2|   553     255     51M     *       0       0       AAAATTTCCGTGCGGGCCGTGGGGGGCTGGCGGCAACTGGGGGGCCGCAGA     BBBFFFFFFFFFFIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFF   XA:i:0  MD:Z:51 NM:i:0
    # QNAME FLAG RNAME POS(leftmost) MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    misCount = 0
    lineStr = line
    line = line.split()
    if len(line) >= 13:
        m = line[12].split(':')
        # if(line[1]=='0' or line[1]=='16') and len(m[2]) > 2: #mismatch detected
        # if (line[1]=='0' or line[1]=='16'):
        # if m[1] != 'Z':

        if line[1] == '0':  # postive strand, remove 5' end

            while m[2].startswith('0'):  # have a mismatch
                m[2] = m[2][2:]  # strip first two character
                misCount += 1  # mismatch count +1
            if misCount != 0:
                line[3] = str(int(line[3]) + misCount)  # move 5'start to 3' direction
                line[5] = str(int(line[5].rstrip('M')) - misCount) + 'M'  # change seq length
                line[9] = line[9][misCount:]
                line[10] = line[10][misCount:]
                line[12] = 'MD:Z:' + str(len(line[9]))
            lineStr = '\t'.join(line)
            lineStr = lineStr + '\n'
        elif line[1] == '16':  # negtive strand, remove 3' end
            while (m[2][-1] == '0' and m[2][-2] in "ATCG"):  # have a mismatch
                m[2] = m[2][:-2]  # strip last two character
                misCount += 1  # mismatch count +1
            if misCount != 0:
                line[5] = str(int(line[5].rstrip('M')) - misCount) + 'M'  # change seq length
                line[9] = line[9][:-misCount]
                line[10] = line[10][:-misCount]
                line[12] = 'MD:Z:' + str(len(line[9]))
            lineStr = '\t'.join(line)
            lineStr = lineStr + '\n'
    return (lineStr, misCount)