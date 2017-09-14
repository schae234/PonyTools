import pytest

from ponytools.AxiomCalls import AxiomCalls
from ponytools.AxiomAnnot import AxiomAnnot
from ponytools.VCF import VCF

import pandas as pd

from locuspocus import Fasta
from locuspocus.Fasta import Chromosome


annot_text = '''"Probe Set ID","Affy SNP ID","dbSNP RS ID","Chromosome","Physical Position","Strand","Flank","Allele A","Allele B","cust_chr","cust_id","cust_pos","ChrXpseudo-autosomal region 1","ChrX pseudo-autosomal region 2","Genetic Map"
"probe1","SNP1","---","1","100","+","TAATCCAAGCCGGCAAAGGCGAAGCTCTGAGGATC[C/G]GCGCCATCCTGCGCTCCCTCGTGCCCCTGGACGAC","A","C","chr1","MNEc_2_1_100_VIP","100","---","---","---"
"probe2","SNP2","---","2","2000","+","GCTTAAAGCAGGCGAGGCCCAGGGAGAGGATAGAA[C/T]GGGTCCTGGCAGCATGACACACGGCCTTGTAACGT","C","G","chr2","MNEc_2_2_2000_VIP","2000","---","---","---"
"probe3","SNP3","---","3","30000","+","CCATAACCATCACCATGATAGGTGTGGTTCTCTTT[A/G]ATTTTGCTGCTGACTTCATTGATGGGCCCATCAAA","G","T","chr3","MNEc_2_3_30000_VIP","30000","---","---","---"
"probe4","SNP4","---","4","400000","+","CAGCCGCTCCTCGGCCTCGGCCGCCGACCGAACTT[A/C]GGCAGAGCCCGAGAGCCTCGTGTTGCCCTCCAACG","T","A","chr4","MNEc_2_4_400000_VIP","400000","---","---","---"'''

@pytest.fixture()
def smpl_annot(tmpdir):
    tmpfile = tmpdir.join('smpl_annot.csv')
    tmpfile.write(annot_text)
    annot = AxiomAnnot.from_file(tmpfile.strpath)
    return annot

@pytest.fixture()
def smpl_fasta():
    ''' A simple fasta that agrees with smpl_annot'''
    fasta = Fasta()
    chr1 = Chromosome('A'*500000)
    chr2 = Chromosome('C'*500000)
    chr3 = Chromosome('G'*500000)
    chr4 = Chromosome('T'*500000)
    fasta.add_chrom('chr1',chr1)   
    fasta.add_chrom('chr2',chr2)   
    fasta.add_chrom('chr3',chr3)   
    fasta.add_chrom('chr4',chr4)   
    return fasta

@pytest.fixture()
def alt_fasta():
    ''' A simple fasta that has all alt calls for smpl_annot'''
    fasta = Fasta()
    chr1 = Chromosome('C'*500000)
    chr2 = Chromosome('G'*500000)
    chr3 = Chromosome('T'*500000)
    chr4 = Chromosome('A'*500000)
    fasta.add_chrom('chr1',chr1)   
    fasta.add_chrom('chr2',chr2)   
    fasta.add_chrom('chr3',chr3)   
    fasta.add_chrom('chr4',chr4)   
    return fasta

@pytest.fixture()
def smpl_call():
    df = pd.DataFrame(
        {
            'Sample1':[0,0,0,2],
            'Sample2':[0,0,0,2],
            'Sample3':[0,0,1,2],
            'Sample4':[0,1,1,2],
            'Sample5':[0,1,1,2]
        },
        index=['probe1','probe2','probe3','probe4']
    )
    return AxiomCalls(df)

vcf_text = '''##fileformat=VCFv4.0
##source=PonyTools:v0.1.2
##AxiomCallFiles=
##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##INFO=<ID=MNEc,Number=0,Type=Flag,Description="Variant is a MNEc SNP">
##INFO=<ID=axiomprobe,Number=1,Type=String,Description="Axiom Probe id for variant">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3 Sample4 Sample5
chr1    100     MNEc_2_1_100_VIP        A       C       .       PASS    MNEc=1;axiomprobe=probe1        GT:GQ   0/0:42  0/0:42  0/0:42  0/0:42  0/0:42
chr2    2000    MNEc_2_2_2000_VIP       C       G       .       PASS    MNEc=1;axiomprobe=probe2        GT:GQ   0/0:42  0/0:42  0/0:42  1/0:42  1/0:42
chr3    30000   MNEc_2_3_30000_VIP      G       T       .       PASS    MNEc=1;axiomprobe=probe3        GT:GQ   0/0:42  0/0:42  1/0:42  1/0:42  1/0:42
chr4    400000  MNEc_2_4_400000_VIP     T       A       .       PASS    MNEc=1;axiomprobe=probe4        GT:GQ   1/1:42  1/1:42  1/1:42  1/1:42  1/1:42'''


real_vcf_text = '''##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3'''

@pytest.fixture
def smpl_vcf(tmpdir):
    tmpfile = tmpdir.join('smpl.vcf')
    tmpfile.write(vcf_text)
    vcf = VCF(tmpfile.strpath,force=True)
    return vcf


@pytest.fixture
def test_vcf(tmpdir):
    tmpfile = tmpdir.join('test.vcf')
    tmpfile.write(real_vcf_text)
    vcf = VCF(tmpfile.strpath,force=True)
    return vcf
