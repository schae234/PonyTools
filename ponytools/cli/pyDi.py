#!/usr/bin/env python3.4

import argparse
import os
import sys

import ipdb as pdb
import pandas as pd
import numpy as np
import logging

from ponytools.VCF import VCF
from ponytools.Variant import Fst
from collections import defaultdict
from itertools import combinations as comb



def pyDi(args):
    
    #Generate a VCF object
    vcf = VCF(args.vcf)

    log = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
                    '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    if args.debug is True:
        log.setLevel(logging.DEBUG)

    
    ind_map = defaultdict(list)
    # Generate a dictionary of ids
    for ind_file in args.inds:
        # Extract basename
        base = os.path.basename(ind_file)
        group_name = os.path.splitext(base)[0]
        # Read in ind files line by line and create dict
        with open(ind_file,'r') as IN:
            for line in IN:
                ind_id = line.strip()
                ind_map[group_name].append(vcf.sample2index(ind_id))
    # grab pairwise combinations of groups
    group_combs = list(comb(ind_map.keys(),2))

    # Allocate a massive SNP x group_comb table
    fst_table = list()
    snp_table = list()
    for i,snp in enumerate(vcf.iter_variants()):
        if i % 500000 == 0:
            log.info("Processing SNP {}".format(i))
        if args.debug is True and i == 10000:
                break
        group_fst = list()
        for group_i,group_j in group_combs:
            # Get the group ids for convenience
            ind_i = ind_map[group_i]
            ind_j = ind_map[group_j]
            alt_freq_i = snp.alt_freq(samples_i=ind_i)
            alt_freq_j = snp.alt_freq(samples_i=ind_j)
            if alt_freq_i is not None and alt_freq_j is not None:
                group_fst.append(Fst(alt_freq_i,alt_freq_j))
            else:
                group_fst.append(np.nan)
        fst_table.append(group_fst)
        snp_table.append([snp.id,snp.chrom,snp.pos])
    # Create pandas tables from the snp and Fst tables
    fst_table = pd.DataFrame(fst_table,columns=group_combs)
    snp_table = pd.DataFrame(snp_table,columns=['snp','chrom','pos'])

    # Standard normalize the columns of the Fst table
    log.debug('Raw Fst:\n{}'.format(fst_table))
    log.debug("Column means:\n{}".format(fst_table.mean()))
    log.debug("Column stdev:\n{}".format(fst_table.std()))
    # Normalize over COLUMNS!!!! (axis is 0)
    fst_table = fst_table.apply( lambda col : (col-col.mean())/col.std(), axis=0)
    # Add windows to the snp table
    log.info("Creating windows that are {} bp large".format(args.window_size))
    snp_table['window'] = snp_table.pos // args.window_size
    mean_window_num = snp_table.groupby(['chrom','window']).apply(len).mean()
    std_window_num = snp_table.groupby(['chrom','window']).apply(len).std()
    log.info("Mean number of SNPs per window {}(+/-){}".format(mean_window_num,std_window_num))

    # Filter out windows with less than the minumum number of SNPs
    snp_table = snp_table.groupby(['chrom','window']).filter(lambda df: len(df) >= args.min_snps_per_window)

    for group in ind_map.keys():
        mask = [group in comb for comb in group_combs]
        # Sum the corrected Fst values for combinations in the group
        snp_table[group] = fst_table.ix[:,mask].sum(axis=1)

    def colmean(col):
        if sum(np.logical_not(np.isnan(col))) >= args.min_snps_per_window:
            return np.nanmean(col)
        else:
            return np.nan

    windowed_di = snp_table.ix[:,['chrom','pos','window']+list(ind_map.keys())].groupby(['chrom','window']).agg(colmean)

    log.info("Extreme Di for each breed:")
    for breed in ind_map.keys():
        log.info("\t%s\tMax:%f,\tMin:%f",breed,windowed_di[breed].max(),windowed_di[breed].min())
    
    # output the tsv files
    snp_table.to_csv(args.out+'_'+"_".join(ind_map.keys())+'_SNPWise_Di.tsv',sep='\t')    
    windowed_di.to_csv(args.out+'_'+"_".join(ind_map.keys())+'_windowed_Di.tsv',sep='\t')    

    if args.debug:
        import pdb; pdb.set_trace()


