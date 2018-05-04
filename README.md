PonyTools
=========


[![Build Status](https://travis-ci.org/schae234/PonyTools.svg?branch=master)](https://travis-ci.org/schae234/PonyTools)


Tools for analyzing Variant (VCF) data and the MNEc Affymetrix SNP Chip

Inlcudes python classes for:
+ VCFs
+ Axiom Genotype Calls
+ FASTA files
+ VCF style variants

Includes command line tools to:
+ Convert AxiomCalls to VCFs
+ Conform VCFs for genotype imputation 
+ Calculate Concordance between VCF Files
+ Calculate Precision and Recall using different scores between a test VCF and a reference VCF
+ Plot various summary statistics from VCF files
+ Calculate the Di statistic using VCF files

Installation
------------
**Ponytool was developed using python 3.6**

For best results, install ponytools into a fresh virtual environment.

### Stable Version
The most recent recent stable version is on pypi
```
pip install ponytools
```

### Bleeding Edge version
To install the developmental version from github:
```
# using anaconda virtual env
conda create -n ponytools python=3
source activate ponytools
git clone git@github.com:schae234/PonyTools.git
cd PonyTools
pip install .
```

