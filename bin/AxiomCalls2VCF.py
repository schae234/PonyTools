#!/usr/bin/env python3
import ponytools as pc


import sys
import argparse


def main(args):
    calls = pc.AxiomCalls.from_file(args.calls,sample_file=args.samples)
    annot = pc.AxiomAnnot.from_file(args.annot)
    fasta = pc.Fasta.from_file(args.fasta)
    calls.to_vcf(args.out,annot,Fasta=fasta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert axiom calls to VCF files') 
    parser.add_argument('--calls',help='Axiom Calls File')   
    parser.add_argument('--samples',help='Axiom Sample File')
    parser.add_argument('--annot',help='Axiom Annotations File')   
    parser.add_argument('--fasta',help='Fasta File')   
    parser.add_argument('--out',help='Oupute file name')
    args = parser.parse_args()
    sys.exit(main(args))
