# GRECO motif clustering pipeline

## Requirements
The pipeline is developed  and tested on UNIX-like operating systems, and has the following dependencies:
* *[Python](https://www.python.org/)* (3.X) with the *[Biopython](https://biopython.org/)*, *[genome-tools](https://github.com/jvierstra/genome-tools)*, *[NumPy](https://numpy.org/)*, *[pandas](https://pandas.pydata.org/)*, *[SciPy](https://www.scipy.org/)* and *[tqdm](https://tqdm.github.io/)* libraries
* *[Tomtom](https://meme-suite.org/meme/doc/tomtom.html)* from the *[MEME suite](https://meme-suite.org/meme/index.html)*

## Installation
To install the previous dependencies, I highly recommend using the conda package manager:
```
conda create -n motif-clustering -c conda-forge -c bioconda \
    biopython meme matplotlib numpy pandas scipy seaborn tqdm
```
To install *genome-tools*, however, first activate the conda environment:
```
conda activate motif-clustering
```
Then clone the git repository and install as follows:
```
git clone https://github.com/jvierstra/genome-tools.git
cd genome-tools
python setup.py install
```

## Quick start
The script *motifs2clusters.py* clusters motifs:
```
$ ./motifs2clusters.py -h
usage: motifs2clusters.py [-h] [--out-dir OUT_DIR] motifs_dir

reformats motifs into MEME format (one file per TF).

positional arguments:
  motifs_dir         motifs directory

optional arguments:
  -h, --help         show this help message and exit
  --out-dir OUT_DIR  output directory (default: ./)
  ```
It takes as input a folder containing motifs belonging to one or more TFs and in *.pcm* and/or *.pcm* format. The script outputs one folder per TF. The content of this folder is a series of intermediate files, and a subfolder called *clusters.0.70 *. This folder contains multiple files:
* *cluster-motifs.$NUM.txt* contains the members of cluster *$NUM*
* *cluster-seed.$NUM.txt* contains the seed motif of cluster *$NUM*
* *cluster.$NUM.pdf* provides a visualization of cluster *$NUM* in *.pdf* format
* *cluster.$NUM.png* provides a visualization of cluster *$NUM* in *.png* format
