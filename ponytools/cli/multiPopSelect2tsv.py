#!/usr/bin/env python3

import sys
import numpy as np
from optparse import OptionParser

def main(args):
    parser = OptionParser()
    parser.add_option('--mpopout')    
    parser.add_option('--out',default=None)
    parser.add_option('--min_pops',dest='minpops',default=0,type='int')
    parser.add_option('--max_pops',dest='maxpops',default=float('inf'),type='float')
    parser.add_option('--exclude_required',action='store_true',default=False)
    options,args = parser.parse_args(args)

    if options.out is None:
        options.out = sys.stdout
    else:
        options.out = open(options.out,'w')
    if len(args) > 0:
        parser.error('Incorrect number or arguments: {}'.format(args)) 

    tagfile = args
    seen_tags = set()
    with open(options.mpopout,'r') as IN:
        num_snps = int(IN.readline().strip().split()[-1])
        tags = IN.readline().strip().split()[1:]
        if len(tags) != num_snps:
            raise ValueError('Number of tag snps doesnt match: found {}, alledged: {}'.format(len(tags)),num_snps)
        # we have to do some replace chaining because the authors of multipopselect hate
        # parseable result files
        tags = [x.replace('-','\t').replace('\t','-',1).replace('/','\t') for x in tags]
        for tag in tags:
            id,*pops = tag.split('\t')
            pops = np.array(list(map(int,pops)))
            if sum(pops > 0) >= options.minpops and sum(pops > 0) <= options.maxpops:
                seen_tags.add(id)
                print(tag,file=options.out)
            elif all(pops == 0) and id not in seen_tags and not options.exclude_required:
                print(tag,file=options.out)
        

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
