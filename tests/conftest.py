import pytest

from ponytools.AxiomCalls import AxiomCalls
from ponytools.Fasta import Fasta
from ponytools.AxiomAnnot import AxiomAnnot
from ponytools.Tools import Chromosome
from ponytools.VCF import VCF

import pandas as pd


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

@pytest.fixture
def smpl_vcf(tmpdir):
    tmpfile = tmpdir.join('smpl.vcf')
    tmpfile.write(vcf_text)
    vcf = VCF(tmpfile.strpath,force=True)
    return vcf
