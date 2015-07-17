'''
        Tests
''' 

import pytest
import ponytools as pc

@pytest.fixture
def smpl_var():
    return pc.Variant(
        'chr1',100,'SNP1',
        'A','C',
        qual=100,fltr='PASS',
        info='.',fmt='GT:GQ',
        genos=['0/0:42','0/1:42','1/1:42']
    )

def test_variant_init(smpl_var):
    var = smpl_var
    assert var.chrom == 'chr1'
    assert var.pos == 100
    assert var.id == 'SNP1'
    assert var.ref == 'A'
    assert var.alt == 'C'
    assert var.qual == 100
    assert var.fltr == 'PASS'
    assert var.info == '.'
    assert var.fmt == 'GT:GQ'
    assert len(var.genos) == 3

def test_variant_from_str(smpl_var):
    str_var = pc.Variant.from_str(
        'chr1\t100\tSNP1\tA\tC\t100\tPASS\t.\tGT:GQ\t0/0:42\t0/1:42\t1/1:42'
    )
    assert str_var == smpl_var

def test_init_from_dict(smpl_var):
    dict_var = pc.Variant.from_dict(
        {   'chrom' : 'chr1',
            'pos'   : 100,
            'id'    : 'SNP1',
            'ref'   : 'A',
            'alt'   : 'C',
            'qual'  : 100,
            'fltr'  : 'PASS',
            'info'  : '.',
            'fmt'   : 'GT:GQ',
            'genos' : ['0/0:42','0/1:42','1/1:42']
        }
    )
    assert dict_var == smpl_var
    
def test_variant_eq(smpl_var):
    assert smpl_var == smpl_var

def test_conform(smpl_var):
    var = smpl_var
    var = var.conform('C')
    assert var.ref == 'C'
    assert var.alt == 'A'
    assert var.genos[0] == '1/1:42'
    assert var.genos[1] == '1/0:42'
    assert var.genos[2] == '0/0:42'

def test_already_conformed(smpl_var):
    var = smpl_var
    var = var.conform('A')
    assert var.ref == 'A'
    assert var.alt == 'C'
    assert var.genos[0] == '0/0:42'
    assert var.genos[1] == '0/1:42'
    assert var.genos[2] == '1/1:42'

def test_add_info(smpl_var):
    smpl_var.add_info('foo','bar')
    assert 'foo=bar' in smpl_var.info.split(';')

