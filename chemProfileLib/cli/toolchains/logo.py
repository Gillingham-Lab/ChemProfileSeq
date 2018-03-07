import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argh
from typing import Union

from chemProfileLib import Configuration
from chemProfileLib.analytics import filter, shannon
from chemProfileLib.helpers import plots

@argh.arg("--model", help="Either the subtraction model or the common peaks",
          choices=["subtraction", "pure", "reactivities"])
@argh.arg("--base0", help="Logo if the modification happened at base0",
          choices=["A", "C", "G", "U"])
@argh.arg("--basep1", help="Logo if the modification happened at base +1",
          choices=["A", "C", "G", "U"])
@argh.arg("--basen1", help="Logo if the modification happened at base -1",
          choices=["A", "C", "G", "U"])
@argh.arg("--logoOf", help="Logo of what kind of modification?",
          choices=["all", "stops", "dels", "ins", "muts"])
def logo(
        configuration: "Path to the configuration file.",
        sample_id: "Sample ID, or a comma separated list",
        genes: "List of genes to calculate logo" = None,
        model = "subtraction",
        base0 = None,
        basep1 = None,
        basen1 = None,
        around = 3,
        logoOf = "stops"
):
    config = Configuration(configuration)
    base0 = "T" if base0 == "U" else base0
    basep1 = "T" if basep1 == "U" else basep1
    basen1 = "T" if basen1 == "U" else basen1

    infile = config.extend_output_path("{sample}.Modcount.data")

    outfile = os.path.join(config.extend_output_path("images"), "Logo.{experiment}.{threshold}.png")

    try:
        control, samples = [s.strip() for s in sample_id.split(":")]
        control = [s for s in control.split(",")]
        treated = [s for s in samples.split(",")]
        samples = control + treated
    except ValueError:
        if model == "pure":
            samples = [sample_id]
        else:
            print("You must give at least 1 control and 1 sample if you are not using the pure model.")
            exit(1)

    features = {}
    first_sample_done = False
    samples_done = 0
    features_skipped_lowpass = 0

    logo_method = {
        "subtraction": logo_subtraction,
        "pure": logo_pure,
        "reactivities": logo_reactivities,
    }

    for sample in samples:
        print("Loading sample {}".format(sample))
        samples_done += 1
        with open(infile.format(sample=sample), "rb") as fh:
            while True:
                try:
                    gene_info, gene_data = pickle.load(fh)
                    gene_id = gene_info["gene_id"]
                except EOFError:
                    first_sample_done = True
                    break

                if first_sample_done is False:
                    features[gene_id] = Gen(gene_info)

                if gene_id in features:
                    if gene_data[:,1].mean() >= 1:
                        #print(gene_data[:,1].sum(), len(gene_data[:,1]))
                        try:
                            features[gene_id].append_data(sample, gene_data)
                        except ValueError:
                            continue
                    else:
                        features_skipped_lowpass += 1

        # Clean up
        feature_ids = list(features.keys())
        for id in feature_ids:
            if len(features[id]) < samples_done:
                del features[id]

        print("\tCurrent number of gene Features: {}".format(len(features)))
        print("\tSkipped because of lowpass filter: {}".format(features_skipped_lowpass))
        features_skipped_lowpass = 0

    # Motifs
    motifs = Motifs()

    feature_ids = list(features.keys())
    for id in feature_ids:
        feature = features[id]
        if model == "pure":
            motifs.add(feature.apply(logo_method[model], samples, samples, around=around, base0=base0, basep1=basep1, basen1=basen1))
        else:
            motifs.add(feature.apply(logo_method[model], control, treated, around=around, base0=base0, basep1=basep1, basen1=basen1))

    #motifs.display()
    entropy = motifs.get_entropy()
    entropy_length = None
    entropy_middle = None

    for key in entropy:
        entropy_length = entropy[key].shape[0] if entropy_length is None else entropy_length
        entropy_middle = entropy[key].shape[0]//2 if entropy_middle is None else entropy_middle

        if logoOf != "all" and logoOf == key:
            print("{what}\tA\tC\tG\tU".format(what=key))

            for i in range(0, entropy_length):
                formats = [entropy[key][i,x] for x in range(4)]
                formats = [i-entropy_middle] + formats
                print("{}\t{:.03f}\t{:.03f}\t{:.03f}\t{:.03f}".format(*formats))

            print(".\n")


bases = "ACGT"
base2index = {"A":0,"C":1,"G":2,"T":3}

class Gen():
    data = {}
    length = None

    def __init__(self, gene_info):
        self.gene_info = gene_info
        self.data = {}

    def append_data(self, id, data):
        if self.length is None:
            self.length = len(data)
        elif self.length != len(data):
            raise ValueError("Data must have the same length.")

        self.data[id] = data[:,0:5]
        self.generate_sequence(data[:,4:8])

        # mutational analysis
        seq_data = data[:,4:8]
        for i in range(0, 4):
            seq_data[seq_data[:, i] > data[:, 0] * 0.1, i] = 0

        self.data[id][:,4] = seq_data.sum(axis=1)

    def generate_sequence(self, data):
        sequence = ""
        sequence_data = data.argmax(axis=1)
        pos = 0

        for index in sequence_data:
            base = bases[index]
            sequence += base
            pos += 1

        self.sequence = sequence

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        data = None
        count = 0
        if isinstance(item, list):
            for i in item:
                if i in self.data:
                    if data is None:
                        data = self.data[i]
                        count+=1
                    else:
                        data += self.data[i]
                        count+=1
                else:
                    raise Exception("{} is not found in this gene".format(item))
            data = data / count
        else:
            if item in self.data:
                data = self.data[item]
            else:
                raise Exception("{} is not found in this gene".format(item))

        return data

    def get_motifs(self, countMatrix, around: int = 3, base0: str = None, basep1: str = None, basen1: str = None) -> np.ndarray:
        assert len(countMatrix) == len(self.sequence)
        motif_length = around * 2 + 1

        # Get motifs
        motifs = {}
        for i in range(len(countMatrix)):
            if (i < around) or i > (len(countMatrix)-around):
                continue

            motif = self.sequence[i-around:i+around+1]
            center_base = self.sequence[i]

            # If base0 is given, only look at sequence logos containing that letter at pos 0.
            if base0 is not None and center_base != base0:
                continue
            if basep1 is not None and self.sequence[i+1] != basep1:
                continue
            if basen1 is not None and self.sequence[i-1] != basen1:
                continue

            if countMatrix[i] > 0:
                if motif not in motifs:
                    motifs[motif] = 0

                motifs[motif] += countMatrix[i]

        # Calculate base distribution
        distribution = np.zeros([motif_length, 4])
        for motif in motifs:
            count = motifs[motif]
            j = 0
            for base in motif:
                distribution[j,base2index[base]] += count
                j+=1

        return distribution

    def apply(self, logomethod, control, treated, around: int = 3, base0: str = None, basep1: str = None, basen1: str = None) -> dict:
        return logomethod(self, control, treated, around, base0, basep1, basen1)

class Motifs():
    motifs = {}

    def add(self, motifs: dict):
        for key in motifs:
            if key in self.motifs:
                self.motifs[key] += motifs[key]
            else:
                self.motifs[key] = motifs[key]

    def display(self):
        print(self.motifs)

    def get_entropy(self) -> dict:
        ret = {}
        for key in self.motifs:
            ret[key] = shannon.calc_entropy(self.motifs[key])

        return ret


def logo_common(gene, control, treated, around: int = 3, base0: str = None, basep1: str = None, basen1: str = None) -> dict:
    raw_control = gene[control][:,1:5]
    raw_treated = gene[treated][:,1:5]

    control = filter.percentile_threshold(raw_control, percentile=60)
    treated = filter.percentile_threshold(raw_treated, percentile=60)

    result = filter.remove_common(treated, control)
    result = filter.weigh(result, weighs=[1.0, 0.0, 0.0])

    s = gene.get_motifs(result, around=around, base0=base0, basep1=basep1, basen1=basen1)

    return {"common": s}

def logo_subtraction(
        gene,
        control: Union[str, list, set],
        treated: Union[str, list, set],
        around: int = 3,
        base0: str = None,
        basep1: str = None,
        basen1: str = None) -> dict:
    raw_control = filter.normalize(gene[control])
    raw_treated = filter.normalize(gene[treated])

    corrected = filter.subtract(raw_treated, raw_control)
    corrected = filter.weigh_for_average_reads(corrected, raw_treated)

    return {
        "stops": gene.get_motifs(corrected[:, 1], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "dels": gene.get_motifs(corrected[:, 2], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "ins": gene.get_motifs(corrected[:, 3], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "muts": gene.get_motifs(corrected[:, 4], around=around, base0=base0, basep1=basep1, basen1=basen1)
    }

def logo_reactivities(
    gene,
    control: Union[str, list, set],
    treated: Union[str, list, set],
    around: int = 3,
    base0: str = None,
    basep1: str = None,
    basen1: str = None) -> dict:

    raw_control = filter.normalize(gene[control])
    raw_treated = filter.normalize(gene[treated])

    subtracted = filter.subtract(raw_treated, raw_control)

    subtracted[subtracted < np.percentile(subtracted, 8)] = 0
    #subtracted = filter.normalize_special(subtracted, 98, 92)
    #subtracted[subtracted < 0.6] = 0

    subtracted = filter.weigh_for_average_reads(subtracted, raw_treated)

    return {
        "stops": gene.get_motifs(subtracted[:, 1], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "dels": gene.get_motifs(subtracted[:, 2], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "ins": gene.get_motifs(subtracted[:, 3], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "muts": gene.get_motifs(subtracted[:, 4], around=around, base0=base0, basep1=basep1, basen1=basen1)
    }


def logo_pure(
        gene,
        control: Union[str, list, set],
        treated: Union[str, list, set],
        around: int = 3,
        base0: str = None,
        basep1: str = None,
        basen1: str = None) -> dict:
    treated = filter.normalize(gene[treated])

    return {
        "stops": gene.get_motifs(treated[:, 1], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "dels": gene.get_motifs(treated[:, 2], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "ins": gene.get_motifs(treated[:, 3], around=around, base0=base0, basep1=basep1, basen1=basen1),
        "muts": gene.get_motifs(treated[:, 4], around=around, base0=base0, basep1=basep1, basen1=basen1)
    }


def filter_common(
        genes,
        index_control: Union[str, list, set],
        index_sample: Union[str, list, set]
):
    stp_motifs = None

    for gene_id in genes:
        gene = genes[gene_id]

        raw_control = gene[index_control][:,1:5]
        raw_sample = gene[index_sample][:,1:5]

        control = filter.percentile_threshold(raw_control, percentile=60)
        sample = filter.percentile_threshold(raw_sample, percentile=60)

        result = filter.remove_common(sample, control)
        result = filter.weigh(result, weighs=[1.0,1.0,1.0])

        s = gene.get_motifs(result)

        stp_motifs = s if stp_motifs is None else stp_motifs + s

    print(stp_motifs)

    return stp_motifs


def filter_subtraction(
        genes,
        index_control: Union[str, list, set],
        index_sample: Union[str, list, set]
):
    stp_motifs = None
    del_motifs = None
    ins_motifs = None
    mut_motifs = None

    for gene_id in genes:
        gene = genes[gene_id]

        try:
            raw_control = gene[index_control]
            raw_sample = gene[index_sample]
        except Exception:
            continue

        for i in range(2, 5):
            raw_control[:,i] = raw_control[:,i] / raw_control[:,0]
            raw_sample[:,i] = raw_sample[:,i] / raw_sample[:,0]


        control = raw_control[:,1:5]
        sample = raw_sample[:,1:5]

        # Normal subtraction
        sample = filter.subtract(filter.normalize(sample), filter.normalize(control))

        # Try to "unnormalize"
        sample[:,0] *= raw_sample[:,1].sum()
        sample[:,1] *= raw_sample[:,2].sum() * raw_sample[:,0]
        sample[:,2] *= raw_sample[:,3].sum() * raw_sample[:,0]
        sample[:,3] *= raw_sample[:,4].sum() * raw_sample[:,0]

        s = gene.get_motifs(sample[:,0])
        d = gene.get_motifs(sample[:,1])
        i = gene.get_motifs(sample[:,2])
        m = gene.get_motifs(sample[:,3])

        stp_motifs = s if stp_motifs is None else stp_motifs + s
        del_motifs = d if del_motifs is None else del_motifs + d
        ins_motifs = i if ins_motifs is None else ins_motifs + i
        mut_motifs = m if mut_motifs is None else mut_motifs + m

    print(stp_motifs)
    print(del_motifs)
    print(mut_motifs)
    print(stp_motifs * del_motifs * mut_motifs)

    return del_motifs


def filter_reactivities(
        genes,
        index_control: Union[str, list, set],
        index_sample: Union[str, list, set]
):
    stp_motifs = None

    for gene_id in genes:
        gene = genes[gene_id]

        raw_control = gene[index_control]
        raw_sample = gene[index_sample]

        """for i in range(2, 5):
            if i == 4:
                continue

            raw_control[:,i] = raw_control[:,i] / raw_control[:,0]
            raw_sample[:,i] = raw_sample[:,i] / raw_sample[:,0]"""

        raw_control[raw_control == np.inf] = 0
        raw_sample[raw_sample == np.inf] = 0


        control = raw_control[:,1:5]
        sample = raw_sample[:,1:5]

        try:
            for i in range(0,4):
                if i == 2:
                    continue

                sample[:,i] = filter.calculate_reactivities(sample[:,i], control[:,i], raw_control[:,0])
        except Exception as e:
            continue

        # Multiply reactivities
        R = sample[:,0] * sample[:,1] * sample[:2]

        # Plot
        plot = plots.StackedPlots("Reactivities of {}".format(gene.gene_info["gene_name"]))
        plot.addTrace(sample[:,0], label="Stops")
        plot.addTrace(sample[:,1], label="Dels")
        plot.addTrace(sample[:,3], label="Muts")
        plot.addTrace(R, label="Multiplikation")
        plot.plot("/home/sauterb/data/Sequencing/BSL02/Output-Test/img/reactivities_{}.png".format(gene.gene_info["gene_name"]))

        s = gene.get_motifs(R * raw_control[:,1])

        stp_motifs = s if stp_motifs is None else stp_motifs + s

        print(stp_motifs)

    print(stp_motifs)

    return stp_motifs


def filter_statistics(
        genes,
        index_control: Union[str, list, set],
        index_sample: Union[str, list, set]
):
    stp_motifs = None
    del_motifs = None
    ins_motifs = None
    mut_motifs = None

    for gene_id in genes:
        gene = genes[gene_id]

        raw_control = gene[index_control][:, 1:5]
        raw_sample = gene[index_sample][:, 1:5]

        im, jm = raw_sample.shape
        p_matrix = np.ones([im, jm])

        for j in range(jm):
            control_max = raw_control[:,j].sum()
            sample_max = raw_sample[:,j].sum()

            for i in range(im):
                test = [[raw_control[i,j], control_max], [raw_sample[i,j], sample_max]]

                if test[0][0] < 5 or test[1][0] < 5:
                    continue

                chi2, p, dof, expected = stats.chi2_contingency(test, correction=True)
                p_matrix[i,j] = p

        raw_sample[p_matrix >= 0.005] = 0
        #use_vals = np.all([0 < p_matrix, p_matrix < 0.005], axis=0)
        #raw_sample[use_vals] *= -np.log2(p_matrix[use_vals])

        # Stops
        s = gene.get_motifs(raw_sample[:,0])
        d = gene.get_motifs(raw_sample[:,1])
        i = gene.get_motifs(raw_sample[:,2])
        m = gene.get_motifs(raw_sample[:,3])

        stp_motifs = s if stp_motifs is None else stp_motifs + s
        del_motifs = d if del_motifs is None else del_motifs + d
        #ins_motifs = i if ins_motifs is None else ins_motifs + i
        mut_motifs = m if mut_motifs is None else mut_motifs + m

    print(stp_motifs)
    print(del_motifs)
    print(mut_motifs)
    print(stp_motifs * del_motifs * mut_motifs)
    return stp_motifs * (del_motifs / del_motifs.max() + 1) * (mut_motifs / mut_motifs.max() + 1)