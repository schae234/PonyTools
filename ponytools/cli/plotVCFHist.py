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

import ponytools as pc

def plotVCFHist(args):
    vcf = pc.VCF(args.vcf)
    if args.sumstat == 'maf':
        alt_freqs = [v.alt_freq() for v in vcf.iter_variants()] 
        alt_freqs = list(filter(None,alt_freqs))
        plt.hist(alt_freqs,bins=args.bins)
        plt.xlabel('Alt Allele Freq.')
        plt.ylabel('Observations')
        if args.title is not None:
            plt.title('{} Alt Allele Dist. (n={})'.format(args.title,len(vcf.samples)))
        if args.out is  None:
            args.out = os.path.split(args.vcf)[1].replace('.vcf','_AltFreq.png')
    elif args.sumstat == 'call_rate':
        call_rate = [v.call_rate() for v in vcf.iter_variants()]
        call_rate = list(filter(None,call_rate))
        plt.hist(call_rate,bins=args.bins)
        plt.xlabel('SNP Call Rate')
        plt.ylabel('Frequency')
        if args.title is not None:
            plt.title('{} Call Rate.'.format(args.title,len(vcf.samples)))
        if args.out is  None:
            args.out = os.path.split(args.vcf)[1].replace('.vcf','_CallRate.png')
    elif args.sumstat == 'inter_distance':
        positions = defaultdict(list)
        for var in vcf.iter_variants():
            alt_freq = var.alt_freq()
            if alt_freq is not None  and alt_freq > 0.01:
                positions[var.chrom].append(var.pos) 
        distances = []
        for chrom,pos in positions.items():
            sorted_pos = sorted(pos)
            for i in range(1,len(pos)):
                distances.append(pos[i]-pos[i-1])
 
        distances = list(filter(None,distances))
        plt.xlabel('Inter SNP Distance (log10)')
        plt.hist(np.log10(distances),bins=args.bins)
        plt.ylabel('Frequency')
        if args.title is not None:
            plt.title('{} InterSNP Distance.'.format(args.title,len(vcf.samples)))
        if args.out is  None:
            args.out = os.path.split(args.vcf)[1].replace('.vcf','_InterDistance.png')
    plt.savefig(args.out)

