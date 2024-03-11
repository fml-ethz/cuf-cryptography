# 🔐 cuf-cryptography - Chemical unclonable functions for decentralized cryptography

- [Overview](#overview)
- [Software Requirements](#software-requirements)
- [Installation Guide](#installation-guide)
- [Data Sources](#data-sources)
- [License](#license)

# Overview
This repository contains the data analysis pipeline, in the form of Jupyter Notebooks and Python files, for the following publication:
> Anne M. Lüscher, Andreas L. Gimpel, Wendelin J. Stark, Reinhard Heckel, Robert N. Grass. Chemical unclonable functions for decentralized cryptography. Manuscript accepted.

Specifically:
- the script `optimization.py` runs a parameter grid search to determine optimal parameters for the Jaccard similarity
- the scripts `downsampling.py` and  `downsampling_constrained.py` perform downsampling to identify the number of reads required for successful discrimination of experiments
- the Python files `fileio.py`, `kmers.py`, `analysis.py`, and `filtering.py` supply methods utilized by the aforementioned scripts
- the Jupyter notebook  `hashing.ipynb` performs analysis using the MinHash and Fuzzy Extractor
- all other Jupyter notebooks are for data curation and visualization


# Software requirements

The data analysis has been tested and performed on Windows 10 using Python 3.9.7. The following Python packages were used and are required: 
```
numpy
pandas
scipy
biopython
pillow
reedsolo
datasketch
plotly
```

# Installation guide
To clone this repository from Github, use
```bash
git clone https://github.com/fml-ethz/cuf-cryptography
cd cuf-cryptography
pip install -r requirements.txt
```
These steps should only take a few seconds to a minute, and no further installation is required. No special hardware or resources are required for use.


# Data sources
Most intermediary files are supplied with this repository, in the form of .csv files. Due to their size, the sequencing data is not supplied with this repository, but is instead deposited on [figshare](https://figshare.com/s/5ace57d8a1c360d45302). To generate filtered sequencing data using BBMap (v38.99), the Python script [in the data directory](/data/runner.py) may be used.



# License
This project is licensed under the GPLv3 license, see [here](LICENSE).
