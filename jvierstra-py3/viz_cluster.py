#!/usr/bin/env python

import argparse
from genome_tools import plotting
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
from pylab import rcParams
import re
matplotlib.use("Agg")
np.seterr(divide="ignore")
rcParams["pdf.fonttype"] = 42
rcParams["svg.fonttype"] = "none"

try:
    from . import get_pwms
except:
    from meme import get_pwms

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    Parses arguments provided through the command line.
    """

    # Initialize
    parser = argparse.ArgumentParser()

    # Arguments
    parser.add_argument("meme_file", metavar="motifs.meme")
    parser.add_argument("cluster_seed", metavar="cluster-seed.txt")
    parser.add_argument("cluster_motifs", metavar="cluster-motifs.txt")

    parser.add_argument(
        "--out-dir",
        default="./",
        help="output directory (default: ./)",
    )
    parser.add_argument(
        "-p", "--prefix",
        help="output files prefix (default: infer from input)",
    )

    return(parser.parse_args())

def viz_cluster(meme_file, cluster_seed, cluster_motifs, prefix=None,
    out_dir="./"):

    pwms = get_pwms(meme_file)

    if prefix is None:
        m = re.search("cluster\-seed\.(\d+)\.txt", cluster_seed)
        prefix = f"cluster.{m.group(1)}"

    cluster_seed = pd.read_csv(cluster_seed, delimiter="\t", header=None)
    cluster_seed.columns = ["cl", "seed_motif", "roffset", "s", "e", "N"]

    cluster_motifs = pd.read_csv(cluster_motifs,  delimiter="\t", header=None)
    cluster_motifs.columns = \
        ["id", "consensus", "strand", "w", "loffset", "roffset", "cluster"]

    N = cluster_motifs.shape[0]
    right = np.max(cluster_motifs["loffset"]+cluster_motifs["w"])

    s = cluster_seed["s"][0]
    e = cluster_seed["e"][0]

    padding = 2

    fig = plt.figure()

    width = (right+padding*2)*0.15
    label_size = 2.0

    fig.set_size_inches(width+label_size, N*0.35)

    gs = gridspec.GridSpec(N, 2, wspace=0, width_ratios=[1, label_size/width])

    plt.subplots_adjust(left=0, top=1, bottom=0)

    for i in range(N):

        ax = fig.add_subplot(gs[i,0])

        m = cluster_motifs.iloc[i,:]

        w = m["w"]
        offset = m["loffset"]
        strand = m["strand"]

        plotting.pwm(pwms[m["id"]].T).render(fig, ax, type="ic", xoffset=offset,
            xlim=(0-padding, right+padding), rc=True if strand=="-" else False)

        ax.axvline(s, color="black", ls='--')
        ax.axvline(e, color="black", ls='--')

        ax.text(0.90, 0.5, m["id"], transform=ax.transAxes, ha="left",
            va="center")

        [ax.spines[loc].set_color("none") for loc in ["top", "right"]]

    plt.savefig(os.path.join(out_dir, f"{prefix}.pdf"))
    plt.savefig(os.path.join(out_dir, f"{prefix}.png"))

    plt.close("all")

def main():

    # Parse arguments
    args = parse_args()

    viz_cluster(args.meme_file, args.cluster_seed, args.cluster_motifs,
        args.prefix, args.out_dir)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()