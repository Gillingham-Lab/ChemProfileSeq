import os
import pickle
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import argh
from typing import Union

from chemProfileLib import Configuration
from chemProfileLib.analytics import filter, shannon
from chemProfileLib.helpers import plots
from chemProfileLib.analytics.stats import cmh, BHcontrol

# Parts of the generation of this enrichment-table is originally from Mod-Seeker
#
#Mod-seq data analysis pipeline
#Copyright (C) 2014  Yizhu Lin
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

def enrichment_table(
    configuration: "Path to the configuration file.",
    sample_id: "Comma separated file identifiers. A colon separates controls from treated, eg DMSO1,DMSO2:DMS1,DMS2",
    genes: "List of genes to calculate logo" = None,
    count_type: "stops or dels" = "stops",
    FDR: "false discovery rate" = 0.05,
    OddsRatio: "Minimum odds ratio" = 1.5
):
    config = Configuration(configuration)


    infile = config.extend_output_path("{sample}.Modcount.data")
    samples = [s.strip() for s in sample_id.split(":")]
    if genes is not None:
        genes = genes.split(",")

    if count_type == "stops":
        outfile = config.extend_output_path("{experiment}.Stops.tsv")
        data_index = 1
    elif count_type == "dels":
        outfile = config.extend_output_path("{experiment}.Dels.tsv")
        data_index = 2
    elif count_type == "muts":
        outfile = config.extend_output_path("{experiment}.Muts.tsv")
        data_index = 4
    else:
        raise Exception("No count_types selected.")

    if len(samples) != 2:
        print("Format MUST be control;treated, c1,c2;t1,t2 or similar.")
        print("Found", samples)
        return

    control = [s.strip() for s in samples[0].split(",")]
    treated = [s.strip() for s in samples[1].split(",")]

    if len(control) != len(treated):
        print("Number of controls must match number of treatments (control {:d}, treated {:d})".format(len(control), len(treated)))
        return
    else:
        print("Using the following samples:")
        print("\tControl: ", *control)
        print("\tTreated: ", *treated)

    must_find = len(control + treated)
    missing = []
    for sample in control + treated:
        filename = infile.format(sample=sample)
        if os.path.exists(filename):
            must_find-=1
            missing += [sample]

    if must_find > 0:
        print("{:d} files have not been found".format(must_find), *missing)
        return

    # Load data
    features = {}
    sample_count = 0
    for sample in control + treated:
        filename = infile.format(sample=sample)

        with open(filename, "rb") as fh:
            while True:
                try:
                    gene_info, gene_data = pickle.load(fh)
                except EOFError:
                    sample_count += 1
                    break

                if genes is not None and gene_info["gene_name"] not in genes:
                    continue

                if gene_info["gene_length"] > 100000:
                    continue

                gene_id = gene_info["gene_id"]
                if sample_count == 0:
                    features[gene_id] = Feature(gene_info)

                if gene_id in features:
                    if sample in control:
                        features[gene_id].add_control(sample, gene_data)
                    else:
                        features[gene_id].add_treated(sample, gene_data)

        if sample_count > 1:
            # Clean up data that is not needed
            ids = list(features.keys())

            for id in ids:
                feature = features[id]

                if len(feature.control) + len(feature.treated) != sample_count:
                    del features[id]
                    continue

    print("Data loaded. Number of genes found is {}".format(len(features)))

    # Get a list of all feature ids so we can modify the dict while running through it.
    ids = list(features.keys())

    # Write data to file
    with open(outfile.format(experiment=sample_id), "w") as fh:
        # Initialise title line
        title = [
            "GeneID",
            "Gene",
            "Position",
            "Nucleotide",
        ]

        # Append the title line with meta information about the used samples
        for c in control:
            title.append("{control} RT {type}".format(control=c, type=count_type))
            title.append("{control} Total RT {type}".format(control=c, type=count_type))
        for t in treated:
            title.append("{treated} RT {type}".format(treated=t, type=count_type))
            title.append("{treated} Total RT {type}".format(treated=t, type=count_type))

        # Append the rest of the title.
        title = title + [
            "Enrichment Ratio (Treated/Control)",
            "CMH p values",
            "p_sig",
            "OR",
            "ORL",
            "ORU",
            "Class",
            "Context"
        ]

        # Write the title
        fh.write("{}\n".format("\t".join(title)))

        # Run through all feature ids
        for id in ids:
            feature = features[id]
            feature.set_zero(32)

            # Skip genes that are too long.
            if len(feature) > 100000:
                print("Skip {}, gene too long".format(feature.feature_info["gene_id"]))
                continue

            # Skip genes that are not everywhere in all controls and treated.
            if len(feature.control) + len(feature.treated) != len(control) + len(treated):
                print("Skip {}".format(feature.feature_info["gene_id"]))
                del features[id]
                continue

            # Report which gene we are running.
            print("Running gen {} (id {}, length {})".format(
                feature.feature_info["gene_name"],
                feature.feature_info["gene_id"],
                len(feature)
            ))

            # Initialise a list of p values and values.
            p_values = []
            values = []

            # Calculate CMH statistics and append the p-value to the list of p_values, and the statistics to the value list
            exceptions = 0
            last_exception = None
            for i in range(len(feature)):
                array = []
                for sample in feature:
                    if i >= sample.shape[0]:
                        print("Feature size mismatch")
                        continue
                    array.append(max(sample[i,data_index], 1))
                    array.append(max(sample[:,data_index].sum(), 1))

                try:
                    chi, p, OR, ORL, ORU = cmh(array)
                    p_values.append(p)
                    values.append((chi, p, OR, ORL, ORU))
                except Exception as e:
                    exceptions += 1
                    last_exception = e

            # Report if anything went wrong before and report the last exception
            if exceptions > 0:
                print("{} exceptions found in statistic calculations.".format(exceptions))
                print("Last exception was", last_exception)
                # Skip gene if anything went wrong
                continue

            # Calculate the significant p_value (it takes the false discovery rate into account)
            p_sig = BHcontrol(p_values, FDR)

            if p_sig is not None:
                for i in range(len(feature)):
                    chi, pval, OR, ORL, ORU = values[i]

                    if pval > p_sig:
                        continue
                    if pval == math.nan:
                        continue

                    # Write output!
                    outline = [
                        feature.feature_info["gene_id"],
                        feature.feature_info["gene_name"],
                        i,
                        feature.sequence[i]
                    ]

                    control_value = 0
                    treated_value = 0
                    sum_value = 0

                    try:
                        for c in list(feature.control.values()):
                            outline = outline + [c[i,data_index], c[:,data_index].sum()]
                            control_value += c[i,data_index]/c[:,data_index].sum()
                            sum_value += c[i,data_index]
                        for t in list(feature.treated.values()):
                            outline = outline + [t[i,data_index], t[:,data_index].sum()]
                            treated_value += t[i,data_index]/t[:,data_index].sum()
                            sum_value += c[i,data_index]
                    except IndexError:
                        continue

                    if sum_value == 0:
                        continue

                    control_value /= len(feature.control)
                    treated_value /= len(feature.treated)
                    ratio = treated_value/control_value

                    if ratio < OddsRatio or ratio == math.nan or ratio == math.inf:
                        continue

                    outline.append(ratio)
                    outline.append(pval)
                    outline.append(p_sig)
                    outline.append(OR)
                    outline.append(ORL)
                    outline.append(ORU)
                    outline.append(5 if math.floor(OR) > 5 else math.floor(OR))
                    outline.append(feature.sequence[i-3:i+3+1])

                    outline = [str(x) for x in outline]
                    fh.write("{}\n".format("\t".join(outline)))

            # Clean up
            del features[id]



bases = "ACGT"
base2index = {"A":0,"C":1,"G":2,"T":3}

def analyse_mutations(data, fi=4, ti=8, into=4, limit=0.05):
    """
    Analyzes base counts and summarizes the 3 bottom as "mutations"
    :param data: The data array
    :param fi: from index, the index of data to be used
    :param ti: to index, the last index of data to be used +1 (ti-fi must be 4)
    :param into: Into which index of data the summarizes mutation should be saved
    :param limit: Applies a threshold above which mutations are considered to be from the second allel
    :return:
    """
    seq_data = data[:,fi:ti]
    for i in range(0, 4):
        seq_data[
            seq_data[:, i] > data[:,0]*limit,
            i
        ] = 0
    data[:,into] = seq_data.sum(axis=1)

    return data


def countUpMutations(sample):
    most = sample[:, 4:8].argmax(axis=1) + 4
    rest = sample[:, 4:8].argsort(axis=1)[:, :-1] + 4
    arange = np.arange(len(most))

    most_values = sample[[[arange], [most]]]

    mutations = sample[[np.transpose([arange, arange, arange]), rest]]
    threshold = mutations < (np.transpose([most_values, most_values, most_values]) * 0.1)

    mutations = mutations.sum(axis=1)

    return mutations


class Feature():
    feature_info = {}
    control = {}
    treated = {}
    sequence = None

    def __init__(self, feature_info):
        self.feature_info = feature_info
        self.control = {}
        self.treated = {}
        self.sequence = None

    def add_control(self, id, data):
        if self.sequence is None:
            self.generate_sequence(data[:,4:8])

        data = analyse_mutations(data)

        self.control[id] = data[:,0:5]

    def add_treated(self, id, data):
        if self.sequence is None:
            self.generate_sequence(data[:, 4:8])

        data = analyse_mutations(data)

        self.treated[id] = data[:, 0:5]

    def set_zero(self, border=32):
        for sample in list(self.control.values()) + list(self.treated.values()):
            sample[0:border] = 0
            sample[-border:] = 0

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
        return len(list(self.control.values())[0])

    def __iter__(self):
        control = list(self.control.values())
        treated = list(self.treated.values())
        ret = []

        for i in range(len(control)):
            ret.append(treated[i])
            ret.append(control[i])

        for r in ret:
            yield r

