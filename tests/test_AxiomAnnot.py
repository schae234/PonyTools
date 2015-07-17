'''
    Tests
'''

import pytest

import ponytools as pc

   

def test_init(smpl_annot):
    assert smpl_annot._annot.shape == (4,11)

