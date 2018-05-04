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
```
# using anaconda virtual env
conda create -n ponytools_env python=3
source activate ponytools_env
```

### Stable Version
The most recent recent stable version is on pypi
```
pip install ponytools
```

### Bleeding Edge version
To install the developmental version from github:
```
pip install -e git+https://github.com/schae234/PonyTools.git#egg=ponytools
```

### Check that installation worked
```
$ ponytools --help

usage: ponytools [-h] [--debug] [--interactive] [--overlook] [--version]
                 Available Commands ...

            ______/``'``'-.                       
           (_   6  \    .^                        
         __ `'.__,  |    `'-.                     
        /_ \  /    /      :`^'                    
      /`/_` \/    /       .'                      
      '/  `'-     |.-'`^. `.                      
    ____                    __              __    
   / __ \____  ____  __  __/ /_____  ____  / /____
  / /_/ / __ \/ __ \/ / / / __/ __ \/ __ \/ / ___/
 / ____/ /_/ / / / / /_/ / /_/ /_/ / /_/ / (__  ) 
/_/    \____/_/ /_/\__, /\__/\____/\____/_/____/  
                  /____/                          

optional arguments:
  -h, --help          show this help message and exit
  --debug             Drop into ipdb when something bad happens.
  --interactive       Initiate an ipdb session right before exiting.
  --overlook          Skip analysis if files produced by --out exists.
  --version           Print version information and exit.

Ponytools. Tools related to analyzing the MNEc SNP Chip:
  Use --help with each command for more info

  Available Commands
    help              Prints this help message
    Axiom2VCF (a2vcf)
                      Convert Axiom calls to a VCF file.
    conformVCF        Conform a VCF based on a FASTA File and MNEc2M
                      Annotations
    sortVCF           Sort a VCF file based on a FASTA file
    compareVCF        compare two VCFs. No need to be sorted or conformed, we
                      will handle that!
    VCFInfo           Print information about a VCF File
    VCFStats          Plot stats/histograms based on VCF values, e.g. alt-
                      freq, call_rate, etc.
    VCFPR             Compare precision vs recall of two VCFs
    Di                Python implementation of DI script
    samples           List and Manipulate Samples
    imputation        Setup the imputation pipeline
    dosageR2          Calculate the Dosage r^2 between reference and imputed
                      SNPs

"You aren't even on the right chromosome!"
```
