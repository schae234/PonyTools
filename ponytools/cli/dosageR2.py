from tqdm import *

def dosageR2(args):
    import ponytools as pc
    import pandas as pd
    import numpy as np
    import os

    from ponytools.Allele import Allele
    from scipy.stats import pearsonr

    import matplotlib.pylab as plt

    referenceVCF = pc.VCF(args.refvcf) 
    imputedVCF = pc.VCF(args.impvcf) 

    if args.out == None:
        args.out = os.path.basename(args.impvcf)

    sampleMask = [referenceVCF.sample2index(x) for x in imputedVCF.samples]
    numRefAlleles = len(referenceVCF.samples) * 2

    refGenerator = (x for x in referenceVCF)
    refBuffer = []

    refDosage = []
    impDosage = []
    refAC = []
    ids = []

    def samePos(x,y):
        if x.chrom == y.chrom and x.pos == y.pos:
            return True
        else:
            return False

    def makeCmp(imp,ref):
        imputed_dosage = [float(x) for x in imp.alleles(format_field='DS')]
        reference_dosage = ref.alleles(samples_i=sampleMask,transform=Allele.vcf2geno)
        AC = int(ref.get_info('AC'))
        if AC > numRefAlleles / 2: #Note: num samples equals half the number of alleles
            AC = numRefAlleles - AC
    
        refDosage.extend(reference_dosage)
        impDosage.extend(imputed_dosage)
        refAC.extend([AC for x in range(len(reference_dosage))])
        ids.extend([ref.id for x in range(len(reference_dosage))])

    for i,impvar in enumerate(tqdm(imputedVCF)):
        # Look in the refBuffer for a matching variant
        found_flag = False
        for j,refvar in enumerate(refBuffer):
            if samePos(impvar,refvar):
                makeCmp(impvar,refvar)
                del refBuffer[j]
                found_flag = True
                break
        if found_flag == True:
            next
        # There is nothing in the refBuffer that matches
        while(True):
            refvar = next(refGenerator)
            if samePos(impvar,refvar):
                makeCmp(impvar,refvar)
                break
            else:
                refBuffer.append(refvar)
            if len(refBuffer) >= 1000:
                print('Couldnt find {} in buffer, searched {}'.format(impvar,len(refBuffer)))
                break
    df = pd.DataFrame({'id':ids,'ref':refDosage,'imp':impDosage,'AC':refAC})
    df.to_csv(args.out+'.dosageR2',index=None)

        
    pccs = df.query('ref >= 0 and AC >= 2').sort_values('AC').groupby('AC').apply(lambda x: pearsonr(x['imp'],x['ref'])[0]).fillna(0)
    plt.semilogx(pccs.index,pccs.values,basex=2)
    plt.xlabel('Minor Allele Count')
    plt.ylabel('R squared')
    plt.title(os.path.basename(args.out))
    plt.savefig(args.out+'.png')
    
