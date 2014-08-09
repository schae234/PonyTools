#!/usr/bin/env python3

import sys

def skipno(filename):
    with open(filename,'r') as IN:
        skip = 0
        while IN.readline().startswith("#"):
            skip += 1
        return skip

def log(message,*formatting):
    print(message.format(*formatting),file=sys.stderr)


def SpecialSNPs():
    special_snps = pd.DataFrame()
    special_path = '/project/mccuelab/rob/MNEc2M/special_snps.txt'
    for i,design in enumerate(['A','B','C']):
        annot = "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Annotation.v1/CsvAnnotationFile_{}.v1.txt".format(design,design,i+1)
        special_snps = pd.concat([special_snps,pd.read_csv(annot,skiprows=skipno(annot),sep=',',quotechar='"')])
    special_snps = special_snps[[x in ['chrX','chrY'] for x in special_snps.cust_chr]]
    special_snps = special_snps[["Affy SNP ID","cust_chr"]]
    special_snps.columns = ['probeset_id','chr']
    special_snps['copy_female'] = [2 if chr == 'chrX' else 0  for chr in special_snps.chr ]
    special_snps['copy_male'] = 1
    special_snps.to_csv(special_path,sep="\t",index=False)
 
