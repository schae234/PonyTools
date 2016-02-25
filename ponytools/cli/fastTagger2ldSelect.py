#! /usr/bin/env python3

import sys
from optparse import OptionParser
import itertools
import numpy as np

def main(args):
    parser = OptionParser()
    parser.add_option('--tagfile',dest='tagfile')
    parser.add_option('--maffile',dest='maffile')
    parser.add_option('--idfile',dest='idfile')
    parser.add_option('--out',dest='out',default=None)

    options,args = parser.parse_args(args)

    if options.out is None:
        options.out = options.tagfile.replace('tagrule.txt','ldselect.txt')
  
    snpids = []
    with open(options.idfile,'r') as IN:
        for i,line in enumerate(IN):
            snpid,pos,a1,a2 = line.strip().split()
            # snp ids are stored by enumeration
            snpids.append(snpid) 

    snpmafs = {}
    with open(options.maffile,'r') as IN:
        for line in IN:
            snpid,pos,a1,a2,maf = line.strip().split()
            snpmafs[snpid] = float(maf)
            

    bins = [] 
    with open(options.tagfile,'r') as IN:
        # islice reads in file two lines at a time
        for l1,l2 in itertools.zip_longest(IN,IN,fillvalue=''):
            k,*lhs = l1.strip().split()
            t,*rhs = l2.strip().split()
            lhs = [snpids[int(x)] for x in lhs]
            rhs = [snpids[int(x)] for x in rhs[::3]]
            # calculate average maf
            avg_maf = np.mean([ snpmafs[snp] for snp in lhs+rhs ])
            if np.isnan(avg_maf):
                avg_maf = 0
            bins.append({
                'total_sites' : int(k)+int(t),
                'average_minor_allele_frequency' : int(avg_maf*100),
                'TagSNPs' : lhs,
                'other_snps' : rhs
            })
    bins.sort(key=lambda x: x['total_sites'],reverse=True)

    with open(options.out,'w') as OUT:
        for i,bin in enumerate(bins):
            print("\n".join([
                'Bin {}\ttotal_sites: {}\taverage_minor_allele_frequency:{}%'.format(i,bin['total_sites'],bin['average_minor_allele_frequency']),
                'Bin {}\tTagSnps: {}'.format(i,' '.join(bin['TagSNPs'])),
                'Bin {}\tother_snps: {}'.format(i,' '.join(bin['other_snps'])),
                ''
            ]),file=OUT)
     

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
