#!/usr/bin/env python3
import sys
from optparse import OptionParser

def main(args):
    parser = OptionParser(usage="usage: %prog --tagfile --idfile")
    parser.add_option('--tagfile')
    parser.add_option('--idfile')
    parser.add_option('--out',default=None)
    options,args = parser.parse_args(args)

    if not options.tagfile:
        parser.error('--tagfile not given')

    if options.out is None:
        options.out = sys.stdout
    else:
        options.out = open(options.out,'w')

    # keep a reference of which array indices are tag SNPs
    tags = set([int(tag.strip()) for tag in open(options.tagfile,'r').readlines()])
    with open(options.idfile,'r') as IN:
        # enumerate id file starting at 0
        for i,snpid in enumerate(IN):
            # if you see a tag SNP print the id
            if i in tags:
                print(snpid.strip(),file=options.out)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:])) 
