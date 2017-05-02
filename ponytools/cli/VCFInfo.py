from ponytools import VCF


def VCFInfo(args):
    vcf = VCF(args.vcf) 
    print(
    '''Information for: {filename}
        Num Samples: {nsamples}
        Num SNPs: {nSNPs}
        Phased: {isPhased}'''.format(
        **{
            'filename' : vcf.vcffile.name,
            'nsamples' : len(vcf.samples),
            'nSNPs' : len(vcf.idmap),
            'isPhased' : vcf.is_phased
        })        
    ) 
