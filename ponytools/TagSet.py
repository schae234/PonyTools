from . import VCF

class TagSet(VCF):
    ''' This class is for a special subset of SNPs within a VCF, namely tSNPs.
        It contains methods specific to tag SNPs and doing fun things with them
        such as '''
    def __init__(self,vcffile,force=False):
        raise NotImplementedError()
        super().__init__(vcffile,force=force)

    def effectiveness(self,genor2_file):
        raise NotImplementedError()
        # read in genor2 file
        gr2 = pd.read_table(genor2_file)
        gr2['POS1_is_Tag'] = [p in self.posmap[c] for c,p in gr2[['CHR','POS1']].itertuples(index=False) ]
        gr2['POS2_is_Tag'] = [p in self.posmap[c] for c,p in gr2[['CHR','POS2']].itertuples(index=False) ]
        # If one is a TAG SNP, it MUST be on the correct chromosome
        len([x for x in set(itertools.chain(*gr2[(gr2['R^2'] >= 0.8) & (gr2.POS1_is_Tag ^ gr2.POS2_is_Tag)][['POS1','POS2']].values))])

    def misclassification(self,imputedVCF,refVCF,call_thresh=0.45):
        raise NotImplementedError()
        mis_frac = list()
        for variant in imputedVCF.iter_variants():
            # skip tag snps
            if variant.pos in self.posmap[variant.chrom]:
                continue
            # Mask out missing values
            ref_var = np.array(refVCF.pos(variant.chrom,variant.pos).genos(transform=Allele.vcf2geno))
            ref_mask = ref_var != -1
            # create a mask on genotypes which are below the call threshold
            imp_var = np.array(variant.alleles(transform=Allele.vcf2geno))
            imp_mask = imp_var > call_thresh
            mask = ref_mask & imp_mask
            agreement = (imp_var[mask] == ref_var[mask])
            if len(agreement) > 0:
                mis_frac.append(
                    sum(agreement)/len(agreement)
                )
        return (np.array(mis_frac).mean(),len(mis_frac),len(refVCF))
 
