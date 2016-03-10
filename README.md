PonyTools
=========

Tools for analyzing Variant (VCF) data and the MNEc Affymetrix SNP Chip

Inlcudes python classes for:
+ VCFs
+ Axiom Genotype Calls
+ FASTA files
+ VCF style variants

Includes command line tools to:
+ Convert AxiomCalls to VCFs
+ Calculate Concordance between VCF Files
+ Calculate Precision and Recall using different scores between a test VCF and a reference VCF
+ Plot various summary statistics from VCF files
+ Calculate the Di statistic using VCF files

Installation
------------
Ponytool was developed using python 3.4

Ponytools is available using the populat python package manager [pip](https://pip.pypa.io/en/stable/installing/)
Using pip, installing ponytools is as easy as:
```
pip install ponytools
```
We recommend installing ponytools in an [virtual environement](https://virtualenv.pypa.io/en/latest/)
*We love [virtualenvwrapper](https://virtualenvwrapper.readthedocs.org/en/latest/)*

Otherwise, requirements for ponytools include:
```
  matplotlib>=1.4.3                                                    
  pandas>=0.16                                                        
  ipython>=4.1.2                                                       
  scipy>=0.17.0                                                     
  ipdb>=0.8 
```
after these packages are installed, you should be able to:
```
$ git clone https://github.com/schae234/PonyTools.git
$ cd PonyTools
$ python3 setup.py install --user
```
