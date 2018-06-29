#!/usr/bin/env python3
import sys
import tempfile
import resource
from optparse import OptionParser
from locuspocus import Fasta

import minus80 as m80


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
    if m80.Tools.available(name=fasta_file,dtype="Fasta"):
        fasta = Fasta(fasta_file)
    else:
        fasta = Fasta.from_file(fasta_file) 
    # Iterate through the chromosome keys and open temp files
    try:
        for chrom in fasta.chrom_names():
            temps[chrom] = tempfile.NamedTemporaryFile('w') 
            chroms.append(chrom)
    except Exception as e:
        log('{}: try increasing the open file limit on your system',e)
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
    log("Outputting sorted chroms")
    if out == None:
        out = vcf_file
        if out.endswith('.vcf'):
            out = out.replace('.vcf','.sorted.vcf')
        else:
            out = out.append('.sorted.vcf')
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
