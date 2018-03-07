import os
import time
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches as patches

from chemProfileLib import Configuration

bases = "ACGT"
colors = {
	"A": "red",
	"C": "cyan",
	"G": "blue",
	"T": "yellow",
}

def plotStops(configuration: "Path to the configuration file.",
        sample_id: "Sample ID",
        gene_id: "Gene ID",
        force_overwrite: "Set to True if you want to overwrite existing files" = False):

    config = Configuration(configuration)
    output_path = config.extend_output_path("images")
    config.force_overwrite = force_overwrite

    if os.path.exists(output_path) == False:
        os.mkdir(output_path)

    out_file = os.path.join(output_path,"{gene}.Contour.{sample}.png")
    modcount_file = "{sample}.Modcount.data"
    modcount_file = config.extend_output_path(modcount_file)

    total_time = time.time()

    samples = sample_id.split(",")

    for sample in samples:
        sample = sample.strip()
        analyze_sample(sample, gene_id, modcount_file, out_file)

    end_time = time.time()

    print("Total running time: {}".format(end_time - total_time))

def analyze_sample(sample, gene_id, infile, outfile):
    infile = infile.format(sample=sample)
    outfile = outfile.format(sample=sample, gene=gene_id)

    with open(infile, "rb") as fh:
        while True:
            try:
                data = pickle.load(fh)
                if data[0]["gene_id"] == gene_id:
                    break
            except EOFError:
                # Reached the end of the file
                data = None
                break

        if data is not None:
            print("{sample}\tGene found".format(sample=sample))
            #print(data[0])
            #print(data[1])

    if data is None:
        print("{sample}\tGene not found.".format(sample=sample))
        return

    gene_data = data[0]
    stop_data = data[1][:,0:8]
    sequ_data = data[1][:,4:8].copy()

    # Prepare mutational analysis
    for i in range(0,4):
        sequ_data[sequ_data[:,i] > stop_data[:,0]*0.1,i] = 0
    muta_data = sequ_data.sum(axis=1)

    stop_data[0:5,:] = 0
    stop_data[-5:,:] = 0
    muta_data[0:5] = 0
    muta_data[-5:] = 0

    fig, axes = plt.subplots(nrows=5, ncols=1, sharex=True)
    axes[0].plot(stop_data[:,0], linewidth=0.5)
    axes[1].plot(stop_data[:,1], linewidth=0.5)
    axes[2].plot(stop_data[:,2], linewidth=0.5)
    axes[3].plot(stop_data[:,3], linewidth=0.5)
    axes[4].plot(muta_data, linewidth=0.5)

    sequence = ""
    sequence_data = stop_data[:, 4:8].argmax(axis=1)
    pos = 0
    for index in sequence_data:
        sequence += bases[index]
        pos+=1

    axes[0].set_title("Stops of {} on gene {}".format(sample, gene_id))
    axes[0].set_ylabel("Coverage")
    axes[1].set_ylabel("Stops")
    axes[2].set_ylabel("Deletions")
    axes[3].set_ylabel("Insertions")
    axes[4].set_ylabel("Mutations")

    for ax in axes:
        ax.get_yaxis().set_label_coords(1.05, 0.5)

    plt.savefig(outfile, dpi=600)
    plt.clf()

    """def get_shannon(sequence, data: np.ndarray, around):
        assert len(sequence) == len(data)
        assert len(data.shape) == 1

        motif_length = around*2+1

        motifs = {}
        # create and count up motifs
        for i in range(len(data)):
            if i < around or i >= len(data)-around:
                continue

            motif = sequence[i-around:i+around+1]
            if motif not in motifs:
                motifs[motif] = 0
            motifs[motif] += data[i]

        # count up base positions
        base_counts = np.zeros([motif_length, 4])
        for motif in motifs:
            for i in range(motif_length):
                base_counts[i,bases.find(motif[i])] += motifs[motif]

        # convert to shanon entropy
        shannon = base_counts
        shannon /= shannon.sum(axis=1)[:, np.newaxis]
        R = 2 + (shannon * np.log2(shannon)).sum(axis=1)
        shannon *= R[:, np.newaxis]
        shannon = np.round(shannon, 5)

        return shannon

    def plot_shannon(shannon):
        f, ax = plt.subplots()
        indices = np.arange(-len(shannon) // 2 + 1, len(shannon) // 2 + 2)

        ax.set_xlabel("Position")
        ax.set_xticks(indices)
        ax.set_xlim(min(indices) - 0.5, max(indices) - 0.5)

        ax.set_ylabel("Bits")
        ax.set_yticks(np.arange(0, 2, 0.1))
        ax.set_ylim((0, max(1, shannon.max())))

        ax.set_title("Shannon entropy")

        for p in range(len(shannon)):
            x_pos = indices[p] - 0.25
            pos_bases = []
            for n in range(len(bases)):
                pos_bases.append((shannon[p, n], bases[n]))
            pos_bases.sort(key=lambda tup: tup[0])

            y_pos = 0
            for base in pos_bases:
                height = base[0]
                ax.add_patch(patches.Rectangle(
                    (x_pos, y_pos),
                    0.5,
                    height,
                    facecolor=colors[base[1]]
                ))
                plt.text(x_pos, y_pos + height, base[1])

                y_pos += height

    shannon = get_shannon(sequence, stop_data[:,1], around)
    plot_shannon(shannon)
    plt.show()"""