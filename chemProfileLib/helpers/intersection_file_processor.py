import multiprocessing as mp
import numpy as np

from .cigar import remove_softclipping


base_index = {
    "A": 4,
    "C": 5,
    "G": 6,
    "T": 7
}

base_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A"
}

temporary_stuff = {}


class Temp:
    count_mutations = False
    count_deletions = False
    count_insertions = False
    count_stops = False
    tmp = {}

    @staticmethod
    def length():
        length = 1
        if Temp.count_mutations:
            length += 4
        if Temp.count_deletions:
            length += 1
        if Temp.count_insertions:
            length += 1
        if Temp.count_stops:
            length += 1

        return length

    @staticmethod
    def i_coverage():
        return 0

    @staticmethod
    def i_stops():
        if Temp.count_stops:
            return 1
        else:
            return None

    @staticmethod
    def i_deletions():
        shift = 0
        if Temp.count_stops:
            shift+=1

        if Temp.count_deletions:
            return 1 + shift
        else:
            return None

    @staticmethod
    def i_insertions():
        shift = 0
        if Temp.count_stops:
            shift += 1
        if Temp.count_deletions:
            shift += 1

        if Temp.count_insertions:
            return 1 + shift
        else:
            return None

    @staticmethod
    def i_bases():
        shift = 0
        if Temp.count_stops:
            shift += 1
        if Temp.count_deletions:
            shift += 1
        if Temp.count_insertions:
            shift += 1

        if Temp.count_mutations:
            return 1 + shift
        else:
            return None


def pool_processor(infile, genome, minReadsPerGene: int, threads: int, countType = None):
    # Prepare read generator
    read_generator = iterateTabFile(infile, genome, minReadsPerGene)

    # prepare pool with general.threads-1 processes in order to save one for saving the data
    p = mp.Pool(max(1, threads))

    for result in p.imap_unordered(process_gene, read_generator):
        yield result

    p.close()

def iterateTabFile(filename, genome, minReadsPerGene):
    if isinstance(filename, str):
        filename = open(filename, "r")

    with filename as fh:
        current_geneid = "_"
        current_line = None
        line_buffer = None
        limit = 9500000

        for line in fh:
            if line.startswith(current_geneid) is False:
                # The line does not start with with the old gene_id
                # We must therefore assume that it is a new one.

                # First, we yield the old stuff (if current_line is not None)
                if line_buffer is not None and len(line_buffer) >= minReadsPerGene:
                    print(
                        len(line_buffer),
                        current_line["gene_id"],
                        current_line["gene_name"],
                        current_line["chr"]
                    ) if len(line_buffer) > limit else None

                    yield  {
                        "gene_id": current_line["gene_id"],
                        "gene_name": current_line["gene_name"],
                        "gene_start": current_line["gene_start"],
                        "gene_stop": current_line["gene_stop"],
                        "strand": current_line["strand"],
                        "type": current_line["type"],
                        "chr": current_line["chr"],
                        "genome_sequence": genome.fetch(
                                current_line["chr"],
                                current_line["gene_start"],
                                current_line["gene_stop"]
                            ).upper(),
                        "reads": line_buffer[:limit],
                    }

                # Then, we collect the new stuff
                current_line = parseLine(line)
                line_buffer = []
                current_geneid = current_line["gene_id"]

            # Append unparsed line
            if current_line["gene_id"] in ["rna69", "gene53975", "gene53976"]:
                continue

            line_buffer.append(line)

def parseLine(line, all=True):
    line = line.split("\t")

    if all:
        columns = {
            "gene_id": line[0],
            "chr": line[1],
            "read_start": int(line[2]),
            "read_end": int(line[3]),
            "qname": line[4],
            "strand": line[5],
            "type": line[6],
            "gene_start": int(line[7]) - 1,
            "gene_stop": int(line[8]),
            "cigar": line[9],
            "read": line[10],
            "gene_name": line[11].strip() if len(line[11].strip()) > 0 else line[0]
        }
    else:
        columns = {
            "read_start": int(line[2]),
            "read_end": int(line[3]),
            "qname": line[4],
            "cigar": line[9],
            "read": line[10],
        }

    return columns

def process_gene(kwargs):
    """
    Processes one gen at a time. Takes in a dictionary called kwargs which matches Gen's constructor.
    :param kwargs: A dictionary matching Gen.__init__
    :return:
    """
    #a = time.time()

    try:
        gen = Gen(**kwargs)
        gen.process_reads()
    except AssertionError:
        return None

    #a = time.time() - a

    """if a > 1:
        print("\t\t\tGene {} with {} reads used {:.02f} seconds ({:.05f} s/read).".format(
            kwargs["gene_id"],
            len(kwargs["reads"]),
            a,
            (a) / len(kwargs["reads"])
        ))"""

    return gen.save()


class Gen():
    gene_id = None
    gene_name = None
    gene_start = 0
    gene_stop = 0
    strand = None
    gene_type = None
    chr = None
    sequence = None
    reads = None
    processed_reads = 0

    def __init__(self, gene_id, gene_name, gene_start, gene_stop, strand, type, chr, genome_sequence, reads):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_start = gene_start
        self.gene_type = type
        self.gene_stop = gene_stop
        self.gene_length = gene_stop - gene_start
        self.chr = chr
        self.strand = strand

        assert len(genome_sequence) == self.gene_length
        self.sequence = genome_sequence

        self.reads = reads
        self.data = np.zeros([self.gene_length, 9], dtype=np.int32)
        self.processed_reads = 0

    def append_read(self, read):
        self.reads.append(read)

    def process_reads(self):
        for read in self.reads:
            new_read = process_read(self.data, read, self.sequence, self.gene_start, self.gene_stop, True if self.strand == "+" else False)
            if new_read is not None:
                self.processed_reads += 1


    def save(self):
        data = self.data if self.strand == "+" else np.flip(self.data, axis=0)

        return {
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "gene_type": self.gene_type,
            "gene_start": self.gene_start,
            "gene_stop": self.gene_stop,
            "gene_length": self.gene_length,
            "chr": self.chr,
            "strand": self.strand,
            "reads": self.processed_reads,
        }, data


def process_read(data: np.ndarray, line: str, sequence: str, gene_start: int, gene_stop: int, pos_strand: bool):
    line = parseLine(line, False)

    read = line["read"]
    cigar = line["cigar"]
    read_start = line["read_start"]
    read_end = line["read_end"]

    new_read, new_cigar = remove_softclipping(read, cigar)

    # Shift chromosome coordinates to transcript coordinates
    start = read_start - gene_start
    end = read_end - gene_start

    # Check sanity
    if start < 0: return
    if end > gene_stop: return
    if end - start > gene_stop - gene_start: return
    if end >= len(sequence): return

    # Do processing
    process_read_cigar(data, read, new_cigar, start, end, pos_strand, sequence[start:end])

    return new_read


def process_read_cigar(data, read, cigar, start, end, positive, sequence):
    # Add Stop
    if positive:
        data[start - 1, 1] += 1
    else:
        data[end, 1] += 1

    pos_sequence = 0
    pos_read = 0
    buf = ""
    numbers = {x: True for x in "0123456789"}

    for char in cigar:
        if char in numbers:
            buf += char
            continue

        length = int(buf)
        buf = ""

        if char == "M":
            # Add coverage data for match
            data[start + pos_sequence: start + pos_sequence + length, 0] += 1
            # Save bases
            print_sequence = ["", ""]
            for read_position in range(length):
                base = read[read_position + pos_read]
                if base == "N":
                    data[start + pos_sequence + read_position, 8] += 1
                else:
                    index = base_index[base]
                    data[start + pos_sequence + read_position, index] += 1

            # (M)atches consume a read and the reference
            pos_sequence += length
            pos_read += length
        elif char == "I":
            # A insertion gets not added to coverage.
            # But it gets added to the insertion line.
            data[max(0, start + pos_sequence):start + pos_sequence + length, 3] += 1

            # (I)nsertions only consume the read, not the reference
            pos_read += length
        elif char == "D":
            # We still put deletions into coverage
            data[max(0, start + pos_sequence): start + pos_sequence + length, 0] += 1

            # And to it's own line, but only the first occurence
            if positive:
                # On positive strand, this is the 3' one
                data[max(0, start + pos_sequence + length - 1): start + pos_sequence + length, 2] += 1
            else:
                # On negative strand, this is the 5' one
                data[max(0, start + pos_sequence): max(0, start + pos_sequence + 1), 2] += 1

            # (D)eletions only consume reference, not the read
            pos_sequence += length
        elif char == "N":
            # Skips only consume reference
            pos_sequence += length
        elif char == "S":
            # Skips only consume read
            pos_read += length