#!/usr/bin/env python3
import sys
import scipy.stats 
from PonyTools.Tools import VCF
from optparse import OptionParser

def main(args):
    parser = OptionParser()
    parser.add_option('--vcf')
    parser.add_option('--num',default=float('inf'),type='int')

    options,args = parser.parse_args(args)

    samples = VCF(options.vcf).samples
    num = min(len(samples),options.num)
    print('\n'.join(scipy.random.choice(samples,num,replace=False)))


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
