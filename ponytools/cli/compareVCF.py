from tqdm import * 

def compareVCF(args):
    import ponytools as pc
    import pandas as pd
    import numpy as np
    import os
    # Create VCF objects
    VCF1 = pc.VCF(args.vcf1)
    VCF2 = pc.VCF(args.vcf2)

    # Get sample lists:
    smpl1 = VCF1.samples
    smpl2 = VCF2.samples
    smpl_both = set(smpl1).intersection(smpl2)
    smask1 = [smpl1.index(s) for s in smpl_both]
    smask2 = [smpl2.index(s) for s in smpl_both]

    var1Buffer = []
    var1Gen = (x for x in VCF1)
    var2Buffer = []

    def samePos(x,y):
        if x.chrom == y.chrom and x.pos == y.pos:
            return True
        else:
            return False

    def makeCmp(var1,var2):
        var1.discordance(var2,smask1,smask2)

    for var2 in tqdm(VCF2,desc="Comparing Variants:"):
        if not var2.is_biallelic:
            next
        found_flag = False
        # Check to see if the variant is in the buffer
        for j,bufvar in enumerate(var1Buffer):
            if samePos(bufvar,var2):
                makeCmp(bufvar,var2)
                del varBuffer[j]
                found_flag = True
                break
        if found_flag == True:
            next
        # There is nothing in the refBuffer that matches
        while(True):
            var1 = next(var1Gen)
            if not var1.is_biallelic:
                next
            if samePos(var1,var2):
                #makeCmp(var1,var2)
                break
            else:
                var1Buffer.append(var1)
            if len(var1Buffer) >= 1000:
                var2Buffer.append(var2)
                break

