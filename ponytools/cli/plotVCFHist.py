#!/usr/bin/env python3.4

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import argparse
import os
import sys

import ipdb as pdb
import pandas as pd
import numpy as np
import logging

from matplotlib import pyplot as plt
from collections import defaultdict
from pysam import VariantFile

import ponytools as pc
import glob

class SampleMask(object):
    def __init__(self, id_file, samples):
        self.name = os.path.basename(id_file).replace('.ids','')
        self.ids = set([line.strip() for line in open(id_file).readlines()])
        self.mask = np.array([x in self.ids for x in samples])
        self.mafs = []

    def process_genos(self,genos):
        # get the correct genos for the mask
        genos = genos[self.mask] 
        self.mafs.append(self.maf(genos))

    @staticmethod
    def maf(genos):
        genos = genos[np.isfinite(genos)]
        if len(genos) == 0:
            return np.nan
        else:
            return sum(genos)/(2*len(genos))

def record_distance(a,b):
    if a[0] != b[0]:
        return np.nan
    else:
        return abs(a[1]-b[1])

def VCFStats(args):
    vcf = VariantFile(args.vcf)
    # get an array of the samples
    samples = np.array(vcf.header.samples)

    # Generate a set of SNPs to analyze
    if args.lst:
        print('Reading in LST file...',end='')
        lst = set(
            tuple(line.strip().split('\t')) \
            for line in open(args.lst,'r')
        )
        print('Done\n')
    
    # generate sample masks
    sample_masks = []
    for id_file in args.ids:
        sample_masks.append(SampleMask(id_file,samples))
    # Initialize some bogus shit 
    prev_pos = ('',0)
    inter_distance = []
    total_maf = []
    ids = []
    chrom = []
    pos = []
    for i,record in enumerate(vcf):
        if args.lst and (record.chrom,str(record.pos)) not in lst:
            continue
        ids.append(record.id)
        # Calculate MAFs
        genos = np.array([
            sum(x.allele_indices) \
            if None not in x.allele_indices \
            else np.nan \
            for x in record.samples.values()
        ])
        for sample in sample_masks:
            sample.process_genos(genos)
        # Calculate Distance
        cur_pos = (record.chrom,record.pos)
        inter_distance.append(
            record_distance(prev_pos,cur_pos)
        )
        chrom.append(record.chrom)
        pos.append(record.pos)
        prev_pos = cur_pos
        # Calculate Total MAF
        total_maf.append(SampleMask.maf(genos))
        # Update log 
        if i % 100000 == 0: 
            print("On {}".format(i))
    # Generate the Data Frame
    data = pd.DataFrame({x.name:x.mafs for x in sample_masks})
    data.insert(0,'ID',ids)
    data.insert(1,'chr',chrom)
    data.insert(2,'pos',pos)
    data.insert(3,'distance',inter_distance)
    data.insert(4,'maf',total_maf)
    # Output Data    
    data.to_csv(args.out+'.tsv',sep='\t')
    return None
