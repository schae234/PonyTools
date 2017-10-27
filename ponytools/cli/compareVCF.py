from tqdm import * 
from ponytools.Exceptions import TriAllelicError,InfoKeyError
from multiprocessing import Pool,Pipe
import os
import time
from datetime import timedelta
import numpy as np

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
    dis = np.zeros(len(xgeno))
    tot = np.zeros(len(xgeno))
    weird = 0
    if xmap != ymap:
        # we need to conform ygenos to xmap
        weird += 1
        #conformed = [xmap.index(y) for y in ymap]
        #ygeno = ['/'.join(conformed[a],conformed[b]) for a,b in map(lambda x: x.split('/'),ygeno)]
    for i,(x,y) in enumerate(zip(xgeno,ygeno)):
        if x.startswith('.') or y.startswith('.'):
            continue
        elif x == y:
            pass
        elif x == y[::-1]:
            pass
        else:
            dis[i] += 1
        tot[i] += 1
    return (dis,tot,weird)

def compare_chroms(tpl):
    chrom,file1,file2,smask1,smask2 = tpl
    var1,seen1pos = parse_chromfile(file1,smask1)
    var2,seen2pos = parse_chromfile(file2,smask2)
    seen_both = seen1pos.intersection(seen2pos)
    vcf1_only = len(seen1pos) - len(seen_both)
    vcf2_only = len(seen2pos) - len(seen_both)
    
    iter1 = (x for x in var1 if x[1] in seen_both)
    iter2 = (x for x in var2 if x[1] in seen_both)

    tot_dis = np.zeros(len(smask1))
    tot_cmp = np.zeros(len(smask1))
    tot_weird = 0
    for x,y in zip(iter1,iter2):
        x_map = [x[3]]+x[4].split(',')
        y_map = [y[3]]+y[4].split(',')
        (dis,cmp,weird) = cmp_geno(x[9:],y[9:],x_map,y_map)

        tot_dis += dis
        tot_cmp += cmp
        tot_weird += weird
     return (tot_dis,tot_cmp,tot_weird,len(seen_both),vcf1_only,vcf2_only)

def compareVCF(args):
    '''
        Compare the discordance between VCFs
    '''
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

    discordant = pd.DataFrame([x[0] for x in chrom_discordances]).sum()
    compared = pd.DataFrame([x[1] for x in chrom_discordances]).sum()
    total_discordance = (discordant.sum()/compared.sum())*100
    weird = sum([x[2] for x in chrom_discordances])
    seen_both = sum([x[3] for x in chrom_discordances])
    seen_1 = sum([x[4] for x in chrom_discordances])
    seen_2 = sum([x[5] for x in chrom_discordances])

    # Print out total stats
    end_time = time.time()
    elapsed = str(timedelta(seconds=(end_time-start_time)))
    print('------------------------')
    print('Elapsed time: {}'.format(elapsed))
    print('Analyzed {} Samples'.format(len(smpl_both)))
    print('Compared {} SNPs'.format(seen_both))
    print('{} SNPs in vcf1 only'.format(seen_1))
    print('{} SNPs in vcf2 only'.format(seen_2))
    print('------------------------')
    print('Discordance: {}%'.format(total_discordance))
    print('Weird: {}'.format(weird))
    print('------------------------')
    if total_discordance > 0:
        print('Offending individuals:')
        ind_discordant = pd.DataFrame(
            [x[0] for x in chrom_discordances],
            columns=smpl_both
        ).sum()
        ind_compared = pd.DataFrame(
            [x[1] for x in chrom_discordances],
            columns=smpl_both
        ).sum()
        ind_percent = (ind_discordant / ind_compared).sort_values(ascending=False) 
        print(ind_percent[ind_percent > 0])

    for v in temps1.values():
        os.remove(v)
    for v in temps2.values():
        os.remove(v)
    return total_discordance,len(smpl_both)
