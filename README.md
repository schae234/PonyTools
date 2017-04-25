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
Ponytool was developed using python 3.4

For an express installation, use the included installtion script. This will create an anaconda
virtual environment
```
git clone git@github.com:schae234/PonyTools.git
cd PonyTools
./install.sh
source activate ponytools
```

If you already have anaconda installed, you can install all the requirements from the conda
`environment.yml` package:
```
conda create  package:
```
git clone git@github.com:schae234/PonyTools.git
cd PonyTools
conda env create -f environment.yml
source activate ponytools
```
