#!/usr/bin/env python

import argparse
import numpy as np
import os
import pandas as pd

from . import get_pwms

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
    parser.add_argument("tomtom_file", metavar="tomtom.txt")
    parser.add_argument("clusters_file", metavar="clusters.txt")
    parser.add_argument(
        "--cluster",
        default=None,
        help="cluster (defaut: all)",
        type=int,
        )
    parser.add_argument(
        "--out-dir",
        default="./",
        help="output directory (default: ./)",
    )

    return(parser.parse_args())

def relative_info_content(pwm):
    p = pwm/np.sum(pwm, axis = 1)[:,np.newaxis]
    ic = 2+np.sum(p*np.nan_to_num(np.log2(p)), axis = 1)
    ric = p*ic[:,np.newaxis]
    return(ric)

def reverse_complement_pwm(pwm):
    baseorder=[3,2,1,0]
    return(pwm[::-1,baseorder])

def pwm_widths(pwms, ids):
    return([pwms[i].shape[0] for i in ids])

def process_clusters(meme_file, tomtom_file, clusters_file, cluster=None,
    out_dir="./"):

    pwms = get_pwms(meme_file)

    sim = pd.read_csv(tomtom_file, header=None, delimiter="\t", skiprows=1,
        comment="#")
    
    clusters = pd.read_csv(clusters_file, header=0, delimiter="\t", index_col=0)

    for cl in sorted(clusters.cluster.unique()):

        if cluster is not None:
            if cluster != cl:
                continue

        cl_motifs = clusters.index[clusters["cluster"] == cl]

        cl_seed_motif = cl_motifs[0]

        i = (sim[0] == cl_seed_motif) & (sim[1].isin(cl_motifs))

        cl_motifs_info = sim.loc[i]

        # motif width
        cl_motifs_info["w"] = cl_motifs_info[8].str.len()

        # cluster left offset
        left = np.min(-cl_motifs_info[2])
        cl_motifs_info["loffset"] = 0 - left - cl_motifs_info[2]

        # cluster right offset
        right = np.max(cl_motifs_info["loffset"] + cl_motifs_info["w"])
        cl_motifs_info["roffset"] = right - cl_motifs_info["w"] \
                                    - cl_motifs_info["loffset"]

        cl_motifs_info["cluster"] = cl

        cl_motifs_info.drop([0, 2, 3, 4, 5, 6, 7], axis=1, inplace=True)
        cl_motifs_info.rename(columns={1: "id", 8: "consensus", 9: "strand"},
                              inplace=True)

        file_name = os.path.join(out_dir, f"cluster-motifs.{cl}.txt")
        cl_motifs_info.to_csv(file_name, sep="\t", header=None, index=False)

        N = cl_motifs_info.shape[0]

        total_ic = np.zeros(right)

        for i in range(N):
            
            m = cl_motifs_info.iloc[i,:]
            
            w=m["w"]
            pwm=pwms[m["id"]]
            offset=m["loffset"]
            strand=m["strand"]

            if strand == "-":
                ic = relative_info_content(pwm)
            else:
                ic = relative_info_content(reverse_complement_pwm(pwm))

            total_ic[offset:offset+w] += np.sum(ic, axis=1)

        cdf = np.cumsum(total_ic) / np.sum(total_ic)
        s = np.where(cdf>0.05)[0][0]
        e = np.where(cdf>0.95)[0][0]+1

        file_name = os.path.join(out_dir, f"cluster-seed.{cl}.txt")
        with open(file_name, "w") as handle:
            handle.write(
                "\t".join(map(str, [cl, cl_seed_motif, right, s, e, N])) + "\n"
            )

def main():

    # Parse arguments
    args = parse_args()

    process_clusters(args.meme_file, args.tomtom_file, args.clusters_file,
        args.cluster, args.out_dir)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()