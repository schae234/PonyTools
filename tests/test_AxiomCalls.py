'''
    Tests
'''
import pytest
import ponytools as pc

def test_init(smpl_call):
    assert smpl_call._calls.shape == (4,5)

def test_to_vcf(tmpdir,smpl_call,smpl_fasta,smpl_annot):
    tmpfile = tmpdir.join('pytest.vcf')
    smpl_call.to_vcf(tmpfile.strpath,smpl_annot,smpl_fasta)
    vcf = pc.VCF(tmpfile.strpath,force=True) 
    assert vcf.pos('chr1',100).alleles() == ['0/0','0/0','0/0','0/0','0/0']

def test_to_conformed_vcf(tmpdir,smpl_call,alt_fasta,smpl_annot):
    tmpfile = tmpdir.join('pytest.vcf')
    smpl_call.to_vcf(tmpfile.strpath,smpl_annot,alt_fasta)
    vcf = pc.VCF(tmpfile.strpath,force=True) 
    assert vcf.pos('chr1',100).alleles() == ['1/1','1/1','1/1','1/1','1/1']
