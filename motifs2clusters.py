#!/usr/bin/env python

import argparse
from Bio import motifs
import os
import pandas as pd
import pathlib
import subprocess as sp
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

from jaspar2others import reformat_motif
utils = __import__("jvierstra-py3")

parser = argparse.ArgumentParser(description="clusters motifs.")
parser.add_argument("motifs_dir", type=pathlib.Path, help="motifs directory")
parser.add_argument("--out-dir", type=pathlib.Path, default="./",
    help="output directory (default: ./)")

args = parser.parse_args()

extensions = set([".pcm", ".ppm"])

tf2paths = {}
for p in args.motifs_dir.rglob("*"):
    if p.is_file():
        if p.suffix in extensions:
            tf = p.name.split("@")[0].split(".")[0]
            tf2paths.setdefault(tf, [])
            tf2paths[tf].append(p)

tf2motifs = {}
for tf in tf2paths:
    tf2motifs.setdefault(tf, [])
    for p in tf2paths[tf]:
        with open(str(p.absolute())) as handle:
            record = motifs.read(handle, "pfm-four-columns")
            record.matrix_id = str(p)
            record.name = tf
            tf2motifs[tf].append(record)

for tf in tqdm(sorted(tf2motifs), total=len(tf2motifs), bar_format=bar_format):
    tf_dir = os.path.join(args.out_dir, tf)
    if not os.path.isdir(tf_dir):
        os.makedirs(tf_dir)
    meme_file = os.path.join(tf_dir, "motifs.meme")
    if not os.path.isfile(meme_file):
        reformat_motif(tf2motifs[tf], "meme", meme_file)
    tomtom_file = os.path.join(tf_dir, "tomtom.txt")
    if not os.path.isfile(tomtom_file):
        cmd = "tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1" +\
            f" {meme_file} {meme_file} > {tomtom_file}"
        sp.run([cmd], shell=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    clusts_file = os.path.join(tf_dir, "clusters.0.70.txt")
    if not os.path.exists(clusts_file):
        utils.hierarchical(tomtom_file, tf_dir)
    clusts_dir = os.path.join(tf_dir, "clusters.0.70")
    if not os.path.isdir(clusts_dir):
        os.makedirs(clusts_dir)
        df = pd.read_csv(clusts_file, header=0, delimiter="\t", index_col=0)
        utils.process_clusters(meme_file, tomtom_file, clusts_file,
            out_dir=clusts_dir)
        for cl in df.cluster.unique():
            clust_seed = os.path.join(clusts_dir, f"cluster-seed.{cl}.txt")
            clust_motifs = os.path.join(clusts_dir, f"cluster-motifs.{cl}.txt")
            utils.viz_cluster(meme_file, clust_seed, clust_motifs,
                out_dir=clusts_dir)