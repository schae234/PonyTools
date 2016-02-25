#!/usr/bin/env python3
import ponytools as pc

import sys

def AxiomCalls2VCF(args):
    calls = pc.AxiomCalls.from_file(args.calls,sample_file=args.samples)
    annot = pc.AxiomAnnot.from_file(args.annot)
    fasta = pc.Fasta.from_file(args.fasta)
    return calls.to_vcf(args.out,annot,Fasta=fasta)
