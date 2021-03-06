'''
    Tests
'''

import pytest
import ponytools as pc

def init_test(smpl_vcf):
    assert smpl_vcf.pos('chr1',100).alleles() == ['0/0','0/0','0/0','0/0','0/0']


def test_pos_index(smpl_vcf):
    assert smpl_vcf.pos('chr1',100).alleles() == ['0/0','0/0','0/0','0/0','0/0']
    assert smpl_vcf.pos('chr2',2000).alleles() == ['0/0','0/0','0/0','1/0','1/0']
    assert smpl_vcf.pos('chr3',30000).alleles() == ['0/0','0/0','1/0','1/0','1/0']
    assert smpl_vcf.pos('chr4',400000).alleles() == ['1/1','1/1','1/1','1/1','1/1']

def test_concord_vcf(smpl_vcf):
    concordance = smpl_vcf.concord_VCF(smpl_vcf)
    assert all([x[1]==1.0 for x in concordance[1:]])

def test_concrd_AxiomCalls(smpl_vcf,smpl_call,smpl_annot):
    concordance = smpl_vcf.concord_AxiomCalls(smpl_call,smpl_annot)
    assert all([x[1]==1.0 for x in concordance[1:]])


