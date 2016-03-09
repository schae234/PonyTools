import ponytools as pc
import pandas as pd
import numpy as np
import matplotlib.pylab as plt

def VCFPR(args):
    test = pc.VCF(args.vcf_test) 
    ref = pc.VCF(args.vcf_ref) 

    # Get a map for the intersection of samples
    sample_map = []
    for sample in set(test.samples).intersection(ref.samples):
        test_index = test.sample2index(sample)
        ref_index = ref.sample2index(sample)
        sample_map.append((sample,test_index,ref_index))
  
    agreement = []
    for ref_var in ref.iter_variants():
        if ref_var in test:
            # extract the other variant
            test_var = test.pos(ref_var.chrom,ref_var.pos)
            # get the rankable score for the test variant
            if args.score.lower() == 'qual':
                score = test_var.qual
            else:
                score = float(test_var.get_info(args.score))
            # Conform the test_var
            test_var.conform(ref_var.ref)
            # Test the samples that overlap
            test_alleles = [test_var.genomap[x] for x in test_var.alleles()]
            ref_alleles =  [ref_var.genomap[x] for x in ref_var.alleles()]
            for sample,test_sample,ref_sample in sample_map:
                if test_alleles[test_sample] == ref_alleles[ref_sample]:
                    agreement.append((sample,test_var.chrom,test_var.pos,score,True))
                else:
                    agreement.append((sample,test_var.chrom,test_var.pos,score,False))
    # Turn results into a dataframe
    agreement = pd.DataFrame(
        agreement,columns=['sample','chrom','pos','score','agreement']
    ).sort_values(by='score', ascending=False)
    # Calculate precision vs recall
    cummulative = []
    precision = []
    for val in agreement.agreement:
        cummulative.append(val)
        precision.append(sum(cummulative)/len(cummulative))
    agreement['precision'] = precision
    # Create the plot
    plt.plot(range(1,len(agreement)+1),agreement.precision)
    plt.ylim(0,1) 
    plt.xlabel('Variant Rank ({})'.format(args.score))
    plt.ylabel('Precision')
    plt.savefig('{}_PR.png'.format(args.out))
    # Output results table
    agreement.to_csv("{}_PR.csv".format(args.out),sep='\t',index=False)
    
