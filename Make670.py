import MNEcTools.Tools as tools
import pandas as pd

def main():
    # Read in the The SNP List
    TSL = pd.read_table("THE_SNP_LIST.tsv",sep=',')
    INFO = pd.DataFrame(tools.read_vcf("/project/mccuelab/DataRepository/Variant/000_SNP/SNP_LIST.0.9.vcf"))
    INFO = INFO.convert_objects(convert_numeric=True)
    TSL = TSL.merge(right=INFO,how='left',left_on=['cust_chr','cust_pos'],right_on=['chrom','pos'])
    
    BEST = TSL[TSL.BestProbeset == 1]
    # get stats for each tissue type and conversion type
    BEST.groupby(['tissue','ConversionType']).apply(len)

    # Split out blood and hair for the best probes:
    BLOOD = TSL[(TSL.BestProbeset == 1)&(TSL.tissue=='blood')]
    HAIR = TSL[(TSL.BestProbeset == 1)&(TSL.tissue=='hair')]
    
    # Compare hair and blood membership
    blood_only = set(BLOOD[BLOOD.ConversionType == 'PolyHighResolution'].snpid).difference(set(HAIR[HAIR.ConversionType == 'PolyHighResolution'].snpid))
    hair_only = set(HAIR[HAIR.ConversionType == 'PolyHighResolution'].snpid).difference(set(BLOOD[BLOOD.ConversionType == 'PolyHighResolution'].snpid))
    
    # Calculate HWE for hair and blood samples seperately
    BLOOD['HWE']=BLOOD[['n_AA','n_AB','n_BB']].apply(lambda x: tools.HW_chi(*x)[1],axis=1)           
    HAIR['HWE']=HAIR[['n_AA','n_AB','n_BB']].apply(lambda x: tools.HW_chi(*x)[1],axis=1)           
    
    # Make final list of:
    #   HighResPoly
    #   HighResMono
    #   NoMinorHomo with HWE > 10e-5
    
    TagSet = pd.concat([
        BLOOD[
            (BLOOD.ConversionType == 'PolyHighResolution') |
            (BLOOD.ConversionType == 'MonoHighResolution') |
            ((BLOOD.ConversionType == 'NoMinorHom') & (BLOOD.HWE > 10e-5))
        ][['probeset_id','cust_chr','cust_pos','cust_id']],
        HAIR[
            (HAIR.ConversionType == 'PolyHighResolution') |
            (HAIR.ConversionType == 'MonoHighResolution') |
            ((HAIR.ConversionType == 'NoMinorHom') & (HAIR.HWE > 10e-5))
        ][['probeset_id','cust_chr','cust_pos','cust_id']]
    ]).drop_duplicates()
    # REad in the Genotype

    blood = AxiomGenos.from_files(
        "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/Samples.csv",
        *["/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/{}/{}/geno/AxiomGT1.calls.txt".format(x,y)
            for x in ['A','B','C'] for y in ['blood']]
    )
    hair = AxiomGenos.from_files(
        "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/Samples.csv",
        *["/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/{}/{}/geno/AxiomGT1.calls.txt".format(x,y)
            for x in ['A','B','C'] for y in ['hair']]
    )
    geno=blood
    geno.add_samples(hair.genos)
    geno.fasta = Fasta.from_file('/project/mccuelab/DataRepository/Fastas/Equus_cab_nucl_wChrUn1_2.fasta')
    annot = annotation.from_files(
        "Axiom_MNEc2M_A_Annotation.r1.csv",
        "Axiom_MNEc2M_B_Annotation.r1.csv",
        "Axiom_MNEc2M_C_Annotation.r1.csv"
    )
    geno.annots = annot
    geno.to_vcf('SNP.vcf')
 
    
    
if __name__ == '__main__':
    sys.exit(main(argv[1:]))
