#!/usr/bin/env python3
import sys
import tempfile
from optparse import OptionParser


def log(message,*formatting):
    print(message.format(*formatting),file=sys.stderr)

def sortVCF(args):
    import sys,os
    vcf_file = args.vcf
    fasta_file = args.fasta
    out = args.out

    headers = list()
    variants = list()
    cur_byte = 0
    chroms = list() 
    temps  = dict()
    log("Sorting {}",vcf_file)
    # Get the chromosome order
    with open(fasta_file,'r') as FASTA:
        for line in FASTA:
            if line.startswith('>'):
                chrom,*info = line.strip().lstrip('>').split()
                log("Found chromosome {}",chrom)
                chroms.append(chrom)
                temps[chrom] = tempfile.NamedTemporaryFile('w') 
    # Get headers and extract positions with file byte offsets
    log("Reading in VCF: {}",vcf_file)
    with open(vcf_file,'r') as VCF:
        for i,line in enumerate(VCF):
            if line.startswith("#"):
                headers.append(line.strip())
            else:
                chrom,pos,*junk = line.split()
                temps[chrom].write(line)
    # close all temp files
    for key,val in temps.items():
        log("flushing tmp file: {}",key)
        val.flush()
    log("soring chroms")
    with open(out,'w') as OUT:
        # print headers
        print("\n".join(headers),file=OUT)
        for chrom in chroms:
            # read in that chroms bullshit
            with open(temps[chrom].name,'r') as CHROM:
                variants =  CHROM.readlines()
                # sort by position
                variants.sort(key=lambda x: int(x.split()[1]))
                log("printing chrom {}",chrom)
                print("".join(variants),file=OUT,end="")
                temps[chrom].close()
