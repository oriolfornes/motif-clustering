#!/usr/bin/env python

import argparse
import numpy as np
import os
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

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
    parser.add_argument("tomtom_file", metavar="tomtom.txt")
    parser.add_argument(
        "--out-dir",
        default="./",
        help="output directory (default: ./)",
    )

    return(parser.parse_args())

def hierarchical(tomtom_file, out_dir="./"):
    """
    https://github.com/jvierstra/motif-clustering/blob/master/hierarchical.py
    """

    # Initialize
    step=0.1
    start=0.5
    end=1.0

    sim = pd.read_csv(tomtom_file, header=None, delimiter="\t", skiprows=1,
        comment="#")

    simsq = sim.pivot(index=0, columns=1, values=4)
    simsq.fillna(100, inplace=True)

    mat = -np.log10(simsq)
    mat[np.isinf(mat)] = 10

    Z = linkage(mat, method="complete", metric="correlation")

    thresholds=np.arange(start, end+step/100, step)

    for thresh in thresholds:
   
        cl = fcluster(Z, thresh, criterion="distance")

        df = pd.DataFrame(cl, index=mat.index, columns=["cluster"])
        df.index.name="motifs"

        df.to_csv(os.path.join(out_dir, "clusters.%0.2f.txt" % (thresh)), sep="\t", header=True)

def main():

    # Parse arguments
    args = parse_args()

    hierarchical(args.tomtom_file, args.out_dir)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()