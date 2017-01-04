from ponytools import VCF
import numpy as np
import math
import sys

def samples(args):
    keep = VCF(args.vcf).samples
    remove = []
    if args.no_permute != True:
        # Then permute
        keep = list(np.random.permutation(keep))
    if args.split != None:
        if args.split > 1:
            args.split = float("0.{}".format(int(args.split)))
        if np.random.randint(2) == 0:
            split_index = math.floor(len(keep)*args.split)
        else:
            split_index = math.ceil(len(keep)*args.split)
        remove = keep[split_index:]
        keep = keep[0:split_index]
    # Print the keep indivs
    if args.out != sys.stdout:
        keep_out = open(args.out+'.keep','w')
        remove_out = open(args.out+'.remove','w')
    else:
        keep_out = args.out
        remove_out = args.out
    print("\n".join(keep),file=keep_out,end='')
    # Print the remove indivs
    print("\n".join(remove),file=remove_out,end='')
