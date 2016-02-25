#!/usr/bin/env python3
import sys
import ponytools as pc
import argparse

from ponytools.Exceptions import TriAllelicError

def main(args):

    vcf = pc.VCF(args.vcf)
    fasta = pc.Fasta.from_file(args.fasta)

    vcf.header['programs'].append('ponytools_conform_vcf')

    with open(args.out,'w') as OUT:
        # print the header
        for key,values in vcf.header.items():
            for val in values:
                print("##{}={}".format(key,val),file=OUT)
        print(vcf.sample_line,file=OUT)
        for i,variant in enumerate(vcf.iter_variants()):
            if i+1 % 10000 == 0:
                print("On SNP {}".format(i),file=sys.stderr)
            # figure out the fasta ref
            fasta_ref = fasta[variant.chrom][variant.pos]
            try:
                variant.conform(fasta_ref)
            except TriAllelicError as e:
                continue
            print(str(variant),file=OUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Python implementation of DI script')
    parser.add_argument('--vcf',action='store',help='VCF file containing all individuals and genotypes')
    parser.add_argument('--fasta',action='store',help='FASTA file containing the reference genotypes')
    parser.add_argument('--out', action='store', default='', type=str, help='Prepend output file names with this.')
    args = parser.parse_args()

    sys.exit(main(args)) 
