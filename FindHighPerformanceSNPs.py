#!/usr/bin/env python3
import pandas as pd
import sys
import os
from optparse import OptionParser

def skipno(filename):
    with open(filename,'r') as IN:
        skip = 0
        while IN.readline().startswith("#"):
            skip += 1
        return skip

def log(message,*vars):
    print("[LOG] ", message.format(*vars),file=sys.stderr)


def main(args):
    parser = OptionParser()
    options,annot_files = parser.parse_args(args)

    for f in annot_files:
        log("Processing {}",os.path.basename(f))
        # read in the annotation file
        annot = pd.read_table(f,quotechar='"',sep=',',skiprows=skipno(f))
        num_probes = annot.groupby('Affy SNP ID').apply(len)
        dual_probe = set(num_probes[num_probes >= 2].index)
        high_performance = annot[[x in dual_probe for x in annot['Affy SNP ID']]]
        # Output probeset IDs
        print("probeset_id")
        for probe in high_performance['Probe Set ID']:
            print(probe) 


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
