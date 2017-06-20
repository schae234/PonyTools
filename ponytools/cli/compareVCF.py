from tqdm import * 
from ponytools.Exceptions import TriAllelicError,InfoKeyError
from multiprocessing import Pool,Pipe
import os
import time
from datetime import timedelta

# print out variants into temp files
def parse_VCF(VCF):
    import tempfile
    temps = dict()
    seen_chroms = set()
    with open(VCF) as IN:
        for line in IN:
            if line.startswith('#'):
                continue
            else:
                chrom,pos,*junk = line.split()
                if chrom not in seen_chroms:
                    temps[chrom] = tempfile.NamedTemporaryFile('w',prefix='ponytools_',delete=False)
                    seen_chroms.add(chrom)   
                temps[chrom].write(line)

    temps = {k:v.name for k,v in temps.items()}
    return (seen_chroms,temps)

def parse_chromfile(filename,sample_mask):
    var = []
    seenpos = set()
    with open(filename,'r') as CHROMFILE:
        for x in CHROMFILE.readlines():
            # chrom,pos,id,ref,alt,qual,fltr,info,fmt,*genotypes 
            x = x.split()
            # if 
            if len(x[3]) > 1 or len(x[4]) > 1 :
                    continue
            gt = x[8].split(':').index('GT')
            var.append(x[0:9]+[x[i].split(':')[gt].replace('|','/') for i in sample_mask])
            seenpos.add(x[1])
    var.sort(key=lambda x:int(x[1]))
    return (var,seenpos)

def cmp_geno(xgeno,ygeno,xmap,ymap):
    dis = 0
    tot = 0
    weird = 0
    if xmap != ymap:
        # we need to conform ygenos to xmap
        weird += 1
        #conformed = [xmap.index(y) for y in ymap]
        #ygeno = ['/'.join(conformed[a],conformed[b]) for a,b in map(lambda x: x.split('/'),ygeno)]
    for x,y in zip(xgeno,ygeno):
        if x.startswith('.') or y.startswith('.'):
            continue
        elif x == y:
            pass
        elif x == y[::-1]:
            pass
        else:
            dis += 1
        tot += 1
    return (dis,tot,weird)

def compare_chroms(tpl):
    chrom,file1,file2,smask1,smask2 = tpl
    var1,seen1pos = parse_chromfile(file1,smask1)
    var2,seen2pos = parse_chromfile(file2,smask2)
    seen_both = seen1pos.intersection(seen2pos)
    
    iter1 = (x for x in var1 if x[1] in seen_both)
    iter2 = (x for x in var2 if x[1] in seen_both)

    tot_dis = 0
    tot_cmp = 0
    tot_weird = 0
    for x,y in zip(iter1,iter2):
        x_map = [x[3]]+x[4].split(',')
        y_map = [y[3]]+y[4].split(',')
        (dis,cmp,weird) = cmp_geno(x[9:],y[9:],x_map,y_map)

        tot_dis += dis
        tot_cmp += cmp
        tot_weird += weird
    return (tot_dis,tot_cmp)

def compareVCF(args):
    import ponytools as pc
    import pandas as pd
    import numpy as np
    import os
    start_time = time.time()
    # Create VCF objects
    VCF1 = pc.VCF(args.vcf1)
    VCF2 = pc.VCF(args.vcf2)
    
    pool = Pool(args.cores)

    # Get sample lists:
    smpl1 = VCF1.samples
    smpl2 = VCF2.samples
    smpl_both = set(smpl1).intersection(smpl2)
    smask1 = [smpl1.index(s)+9 for s in smpl_both]
    smask2 = [smpl2.index(s)+9 for s in smpl_both]
    
    print('Reading in VCFs....')
    ((seen1,temps1),(seen2,temps2)) = pool.map(parse_VCF,[args.vcf1,args.vcf2])

    # get the chroms that are in both
    chroms = seen1.intersection(seen2)
    # create a list of filenames for chroms seen in both
    chroms = [(c,temps1[c],temps2[c],smask1,smask2) for c in chroms]
   
    print('Comparing chromosomes....')
    chrom_discordances = pool.map(compare_chroms,chroms)

    discordant = sum([x[0] for x in chrom_discordances])
    compared = sum([x[1] for x in chrom_discordances])

    # Print out total stats
    print('Discordance: {}'.format(discordant/compared))
    end_time = time.time()
    elapsed = str(timedelta(seconds=(end_time-start_time)))
    print('Elapsed time: {}'.format(elapsed))
