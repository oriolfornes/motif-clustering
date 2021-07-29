#!/usr/bin/env python

import argparse
from genome_tools.plotting import pwm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
from pylab import rcParams
matplotlib.use("Agg")
rcParams["pdf.fonttype"] = 42
rcParams["svg.fonttype"] = "none"

from . import get_pwms

# import sys
# import os.path

# import pandas as pd
# import numpy as np
# import scipy as sp

# import matplotlib
# matplotlib.use('Agg')

# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec

# from pylab import rcParams
# rcParams['pdf.fonttype'] = 42
# rcParams['svg.fonttype'] = 'none'

# from genome_tools.plotting import pwm

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
        prefix = os.path.splitext(os.path.basename(cluster_seed))[0]

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

        #ic=relative_info_content(pwm) if strand=="+" else relative_info_content(reverse_complement_pwm(pwm))

        pwm(pwms[m["id"]].T).render(
            fig, ax, type="ic", xoffset=offset, xlim=(0-padding, right+padding), rc=True if strand=="-" else False
        )

        #stackedbar(ic.T, ax, sorted=True, xoffset=offset)
        ax.axvline(s, color="black", ls='--')
        ax.axvline(e, color="black", ls='--')

        ax.text(0.90, 0.5, m["id"], transform=ax.transAxes, ha="left",
            va="center")

        #ax.xaxis.set_visible(True)

        #ax.set_xlim(0-padding, right+padding)
        #ax.set_ylim(0, 2.25)
        

        
        #ax.set_ylabel("Bits")

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