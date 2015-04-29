#!/usr/bin/env python3 

import argparse
import pandas as pd
import logging
import sys
import numpy as np
from itertools import chain
from collections import Counter,namedtuple,OrderedDict

Ind = namedtuple('Ind',['fam','id','fat','mot','sex','stat'])
MapSNP = namedtuple('MapSNP',['chr','id','cm','pos'])

class smart_dict(dict):
    def __missing__(self, key):
        return key

class Geno(object):
    
    def __init__(self,map,ped):
        self.info = ped.info
        self.map = map
        self.genos = pd.DataFrame(
            ped._genos,
            index = pd.MultiIndex.from_product(
                [[x for x in self.map._snps.keys()],['a1','a2']],names=['snp','allele']
            ),
            columns = self.info.keys() 
        )

    @classmethod
    def from_file(cls,mapfile,pedfile):
        map = Map.from_file(mapfile)
        ped = Ped.from_file(pedfile)
        return cls(map,ped)

class Map(object):

    def __init__(self): 
        self._snps = OrderedDict()
        self._indices = list()

    def add_snp(self,chr,id,cm,pos):
        self._snps[id] = MapSNP(chr,id,cm,pos)
        self._indices.append(id)

    def __iter__(self):
        return (x for x in self._snps.values())

    @classmethod
    def from_file(cls,filename):
        self = cls()
        with open(filename,'r') as IN:
            for line in IN:
                chr,id,cm,pos = line.strip().split()
                self.add_snp(chr,id,cm,pos)
        return self

class Ped(object):
    geno_map = smart_dict({'0':0,'A':1,'C':2,'G':3,'T':4})

    def __init__(self,geno_map=None):
        self.info = OrderedDict()
        self._genos = None
        if geno_map is not None:
            self.geno_map.update(geno_map)

    def __getattr__(self,i):
        return self._genos[i:i+2,:]

    def to_file(self,filename,sep='\t'):
        genos = self._genos.T
        with open(filename,'w') as OUT:
            for i,indiv in enumerate(self.info.values()):
                print(sep.join(map(str,chain(*[indiv,np.array(genos[i]).flatten()]))),file=OUT)

    def to_AB(self):
        abmatrix = np.matrix(np.zeros(self._genos.shape),dtype='U1')
        for i in range(0,len(self._genos),2):
            genos =  self._genos[i:i+2,:]
            # count the alleles
            counts = Counter(np.array(genos).flatten().tolist())
            # get rid of missing counts
            counts.pop(0)
            # Get the AB alleles
            alleles = sorted(
                counts.keys(),
                key = lambda k : counts[k],
                reverse = True
            ) 
            abmatrix[i:i+2,:][genos == alleles[0]] = 'A'
            if len(alleles) > 1:
                abmatrix[i:i+2,:][genos == alleles[1]] = 'B'
        self._genos = abmatrix
        return self

    def add_ind(self,fam,id,fat,mot,sex,stat,genos):
        self.info[id] = Ind(fam,id,fat,mot,sex,stat) 
        genos = [self.geno_map[x] for x in  genos]
        if self._genos is None:
            self._genos = np.matrix(genos).T
        else:
            # genos[...,None] syntax ensures genos is a (...,1) shaped matrix
            genos = np.array(genos)
            self._genos = np.append(self._genos,genos[...,None],axis=1)

    @classmethod
    def from_file(cls,filename):
        self = cls()
        with open(filename,'r') as IN:
            for line in IN:
                fam,id,fat,mot,sex,stat,*genos = line.strip().split()
                self.add_ind(fam,id,fat,mot,sex,stat,genos)
        return self


def main(args):
    parser = argparse.ArgumentParser(description='Process some PED Files.')
    parser.add_argument('--ped', action='append')
    parser.add_argument('--map', action='store')
    args = parser.parse_args(args)

    # Read in map file
    map = Map.from_file(args.map) 
    for pedfile in args.ped:
        print('Working on ped: {}'.format(pedfile))
        ped = Ped.from_file(pedfile)
        ped.to_AB().to_file("AB"+pedfile)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
