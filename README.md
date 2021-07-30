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
The script *motifs2clusters.py* takes as input a folder containing motifs belonging to one or more TFs and in *.pcm* and/or *.pcm* format:
```
$ ./motifs2clusters.py -h
usage: motifs2clusters.py [-h] [--out-dir OUT_DIR] motifs_dir

motif clustering based on Jeff Vierstra's code

positional arguments:
  motifs_dir         motifs directory

optional arguments:
  -h, --help         show this help message and exit
  --out-dir OUT_DIR  output directory (default: ./)
  ```
It returns one folder per TF as output. The contents of such folder are a series of intermediate files (*e.g.* *motifs.meme*, *tomtom.txt*, etc.), and the subfolder *clusters.0.70 * containing multiple files:
* *cluster-motifs.$NUM.txt* contains the members of cluster *$NUM*
* *cluster-seed.$NUM.txt* contains the seed motif of cluster *$NUM*
* *cluster.$NUM.pdf* provides a visualization of cluster *$NUM* in *.pdf* format
* *cluster.$NUM.png* provides a visualization of cluster *$NUM* in *.png* format
