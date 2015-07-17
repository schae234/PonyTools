#!/usr/bin/env python3
import sys
import os
import scipy.stats
from collections import defaultdict
from itertools import chain,groupby
import pandas as pd
import gzip
import pickle
import numpy as np
import multiprocessing 
import time
import math
from subprocess import Popen, PIPE


def skipno(filename,startswith="#"):
    with open(filename,'r') as IN:
        skip = 0
        while IN.readline().startswith(startswith):
            skip += 1
        return skip

def log(message,*formatting):
    print(message.format(*formatting),file=sys.stderr)


class Chromosome(object): 
    '''
    A Chromosome is a lightweight object which maps indices to 
    string positions.
    '''
    def __init__(self,seq): 
        self.seq = str(seq) 
    def __getitem__(self,pos): 
        #chromosomes start at 1, python strings start at 0 
        return self.seq[int(pos)-1]

    def __len__(self):
        return len(self.seq)
 
def HW_chi(n_AA,n_AB,n_BB):
    n = sum([n_AA,n_AB,n_BB])
    p = (2*n_AA+n_AB)/(2*n)
    q = 1-p
    exp_AA = (p**2)*n
    exp_AB = 2*p*q*n
    exp_BB = (q**2)*n
    return scipy.stats.chisquare([n_AA,n_AB,n_BB],[exp_AA,exp_AB,exp_BB],ddof=1)

def VCFFilterByMap(vcf_file,map_file,output_file):
    include = set()
    num_found = 0
    log('reading in map file: {}',map_file)
    with open(map_file,'r') as MAP:
        for line in MAP:
            chrom,id,cm,pos = line.strip().split()
            if not chrom.startswith('chr'):
                chrom = 'chr'+chrom
            include.add((chrom,pos))
    log('reading in vcf file: {}',vcf_file)
    with open(vcf_file,'r') as VCF:
        with open(output_file,'w') as OUT:
            for i,line in enumerate(VCF):
                if i % 100000 == 0:
                    log("on line {}, found {}",i,num_found)
                line = line.strip()
                if line.startswith('#'):
                    print(line,file=OUT)
                else:
                    chrom,pos,*junk = line.split()
                    if (chrom,pos) in include:
                        num_found += 1
                        print(line,file=OUT)

def sortVCF(vcf_file,fasta_file,temp_dir="/tmp",out="sorted.vcf"):
    headers = list()
    variants = list()
    cur_byte = 0
    chroms = list() 
    temps  = dict()
    log("Sorting {}",vcf_file)
    # Get the chromosome order
    with open(fasta_file,'r') as FASTA:
        for line in FASTA:
            if line.startswith('>'):
                chrom,*info = line.strip().lstrip('>').split()
                log("Found chromosome {}",chrom)
                chroms.append(chrom)
                temps[chrom] = open(os.path.join(temp_dir,chrom+'.tmp'),'w')
    # Get headers and extract positions with file byte offsets
    log("Reading in VCF: {}",vcf_file)
    with open(vcf_file,'r') as VCF:
        for i,line in enumerate(VCF):
            line = line.strip()
            if line.startswith("#"):
                headers.append(line)
            else:
                chrom,pos,*junk = line.split()
                print(line,file=temps[chrom])
    # close all temp files
    for key,val in temps.items():
        log("Closing tmp file: {}",key)
        val.close()
    log("soring chroms")
    with open(out,'w') as OUT:
        # print headers
        print("\n".join(headers),file=OUT)
        for chrom in chroms:
            # read in that chroms bullshit
            with open(os.path.join(temp_dir,chrom+'.tmp'),'r') as CHROM:
                variants =  CHROM.readlines()
                # sort by position
                variants.sort(key=lambda x: int(x.split()[1]))
                log("printing chrom {}",chrom)
                print("".join(variants),file=OUT,end="")
                os.remove(os.path.join(temp_dir,chrom+'.tmp'))
            
           
def phasedVCF2FastTagger(phased_VCF,fasta,annot,outprefix=None):
    if not outprefix:
        # use the base name of phased_VCF
        outprefix = os.path.basename(phased_VCF).replace('vcf','fasttagger')
    if 'gz' in phased_VCF:
        vcf_open = gzip.open
        mode = 'rb'
    else:
        vcf_open = open   
        mode = 'r'
    with vcf_open(phased_VCF,mode) as VCF, \
         open(outprefix+'.data_maf','w') as MAP, \
         open(outprefix+'.data_matrix','w') as MAT:
        for line in VCF:
            line = line.decode()
            if line.startswith('#'):
                continue
            chrom,pos,id,ref,alt,qual,filters,info,formats,*genos = line.strip().split()
            info = info.replace('=',';').split(';')
            info = dict(zip(info[0::2],info[1::2]))
            GT_index = formats.split(":").index("GT")
            if ref != fasta[chrom][pos]:
                log('reference bases do not match: {}'," ".join(line.split()[0:9]))
            elif ',' in alt:
                if ',' in id:
                    # Find out if is VIP SNP, skip dual probe SNPs
                    continue
                mafs = [float(x) for x in info['AF'].split(',')]
                # MAFs below 0.005 are probably sequencing errors
                if min(mafs) < 0.005:
                    alt = alt.split(',')[mafs.index(max(mafs))]
                    info['AF'] = max(mafs)
                    genos = [geno.replace('2','.') for geno in genos] # replace the third allele with missing
                else:
                    log("{} is tri-allelic!: {}",id," ".join([chrom,pos,id,ref,alt,qual,filters,";".join(["{}={}".format(k,v) for k,v in info.items()]),formats]))
                    for probe in id.split(','):
                        log("\t{}"," ".join(map(str,annot.table.ix[probe])))
                        log("\t{}",' '.join([geno.split(':')[GT_index] for geno in genos]))
                    continue
            # fastTagger is peculiar about its map files
            if float(info['AF']) > 0.5:
                minor_allele = ref
                major_allele = alt
                maf = 1.0 - float(info['AF'])
            else:
                minor_allele = alt
                major_allele = ref
                maf = float(info['AF'])
            # print line to MAP file
            print("{}\t{}\t{}\t{}\t{}".format(id,pos,major_allele,minor_allele,maf),file=MAP)
            # print corresponding line to matrix file
            print("\t".join(chain(*[geno.split(":")[GT_index].split('|') for geno in genos])),file=MAT)

      

class FastTagger(object):
    binpath = '/project/mccuelab/shared/bin/FastTaggerV2'
    def __init__(self,binpath=None):
        if binpath:
            self.binpath = binpath

    def parfile(self,maf,r2,basename):
        with open(basename+'.par','w') as OUT:
            print("\n".join(
               'data_maf={}.maf'.format(basename),
               'data_matrix={}.matrix'.format(basename),
               'min_maf={}'.format(maf),
               'window_len=100000',
               'min_r2_1={}'.format(r2),
               'min_r2_2={}'.format(r2),
               'min_r2_3={}'.format(r2),
               'max_len=3',
               'max_merge_window_len=100000',
               'min_bin_size=0',
               'max_covered_times=0',
               'mem_size=0',
               'max_tagSNP_num=0.00',
               'model=MMTagger',
               'output=TagSet',
            ),file=OUT)

    def tag(self,vcffile,maf,r2,chrom=None):
        vcf = VCF(vcffile)
        # Create Output Dir
        base = os.path.join("MAF_{}_r2_{}") 
        os.makedirs(base,exist_ok=True)
        # filter down chromosomes
        for chrom in vcf.iter_chroms():
            # make an output for that chrom
            output = os.path.join(base,chrom)
            os.makedirs(output,exist_ok=True)
            if not os.path.exists(os.path.join(output,'{}_{}'.format(chrom,vcf.vcffile.name))):
                p = Popen([
                    
                ])
