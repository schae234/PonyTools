#!/usr/bin/env python3

import os
import pickle
import numpy as np
import gzip
import random

from collections import defaultdict,OrderedDict
from pysam import VariantFile

from .Tools import log
from .Variant import Variant
from .Allele import Allele

import pandas as pd

class VCFHeader(OrderedDict):
    def __init__(self):
        super().__init__()

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            self[key] = list()
            return self[key]

    def __reduce__(self):
        ''' Reduce is needed to pickle '''
        return type(self),tuple()

    def __str__(self):
        return "{}\n{}".format(
           '\n'.join(['##{}={}'.format(key,val) for key,vals in self.items() for val in vals]),
           '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+self.samples)
        )
  
    @classmethod
    def from_file(cls,filename):
        self = cls()
        if filename.endswith('gz'):
            file = gzip.open(filename,'rt')
        else:
            file = open(filename,'r')
        for line in file:
            if line.startswith('##'):
                key,val = line.lstrip('#').rstrip().split('=',1)
                self[key].append(val)
            elif line.startswith("#"):
                samples = line.strip().split()[9:]
                self.samples = samples
            else:
                break
        return self

class VCF(object):
    def __init__(self,vcffile,force=False):
        if vcffile.endswith('gz'):
            self.vcffile = gzip.open(vcffile,'rt')
        else:
            self.vcffile = open(vcffile,'r')
        self.header = VCFHeader.from_file(vcffile)
        # keep track of a bunch of indexes
        self.idmap = {}
        self.posmap = defaultdict(dict)
        self.indexmap = []
        # experimental genotype data frame
        self._genotypes = None
        # load/create indices
        self.index(force=force)

    @property
    def samples(self):
        return self.header.samples

    def __contains__(self,item):
        if isinstance(item,Variant):
            if item.chrom in self.posmap and item.pos in self.posmap[item.chrom]:
                return True
            else:
                return False
        raise KeyError('Cannot test that item')

    @property
    def header_string(self):
       return '\n'.join(['#{}={}'.format(key,val) for key,vals in self.header.items() for val in vals]) 

    @property
    def genotypes(self):
        '''
            This method returns the internal genotype matrix. 

            Notes
            -----
            To save time when initializing the VCF class,
            the matrix is not loaded until the method is
            called the first time, then it is lazy loaded.
        '''
        if self._genotypes is None:
            self._load_genos()
        return self._genotypes

    def _load_genos(self,force=False,transform=Allele.vcf2geno,randomize_missing=False):
        if self._genotypes is not None and force == False:
            return
        log("Loading genotypes for {}",self.vcffile.name)
        alleles = []
        ids = []
        # Use the pysam API
        for i,record in enumerate(VariantFile(self.vcffile.name)): 
            alleles.append(np.array([
                int(sum(x.allele_indices)) \
                if None not in x.allele_indices \
                else np.nan \
                for x in record.samples.values()
            ]))
            ids.append((
                record.chrom,   
                record.pos,
                record.id,
                record.ref,
                record.alts[0]
            ))
            if i % 100000 == 0 and i > 0:
                print("On {}".format(i))
        #for variant in self.iter_variants():  # We only support bi-allelic SNPS right now!!!!
        #    if not variant.is_biallelic:
        #        continue
        #    alleles.append(variant.alleles(transform=transform))
        #    ids.append((variant.chrom,variant.pos,variant.id,variant.ref,variant.alt))
        self._genotypes = pd.DataFrame(
                alleles,
                index=pd.MultiIndex.from_tuples(ids,names=['chrom','pos','id','ref','alt']),
                columns=self.samples
        )

    def index(self,force=False):
        vcf_dirname  = os.path.dirname(self.vcffile.name)
        vcf_basename = os.path.basename(self.vcffile.name)
        index_name = '.' + vcf_basename + '.pdx'
        index_path = os.path.join(vcf_dirname,index_name)
        if (not os.path.exists(index_path) \
                or force \
                or os.path.getmtime(self.vcffile.name) > os.path.getmtime(index_path)
            ):
            # create the index file
            log('Index file does not exist for {}, indexing now'.format(self.vcffile.name))
            cur_byte = 0
            self.vcffile.seek(0)
            for line in self.vcffile:
                if line.startswith('#'):
                    # we need to pass here so that curbyte gets updated below
                    pass
                else:
                    chrom,pos,ids,*junk = line.strip().split()
                    pos = int(pos)
                    for id in ids.split(','): #sometimes there are multiple ids...
                        self.idmap[id] = cur_byte
                    self.posmap[chrom][pos] = cur_byte
                    self.indexmap.append(cur_byte)
                cur_byte += len(line)
            self.vcffile.seek(0)
            # pickle index file for later use
            pickle.dump((self.idmap,self.posmap,self.indexmap),open(index_path,'wb'))
        else:
            # read the index file
            self.idmap,self.posmap,self.indexmap = pickle.load(open(index_path,'rb')) 

    @property
    def sample_line(self):
        return '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+self.samples)

    def iter_variants(self):
        ''' 
            Returns variant generator, for iteration. Should be memory efficient 
        '''
        self.vcffile.seek(0)
        return (Variant.from_str(line) for line in self.vcffile if not line.startswith('#'))

    def rand_variant(self):
        '''
            Returns a random variant from the VCF File
        '''
        rand_index = random.sample(self.indexmap,k=1)[0]
        self.vcffile.seek(rand_index)
        return Variant.from_str(self.vcffile.readline())

    def __iter__(self):
        return self.iter_variants()

    def iter_chroms(self):
        return (chrom for chroms in self.posmap.keys())

    def sample_call_rate(self):
        ''' returns call rate for samples '''
        pass
            
    def sample2index(self,sample):
        ''' return the genotype index for sample '''
        return self.samples.index(sample)

    def __getitem__(self,item):
        if isinstance(item,list):
            try:
                return self.pos(*item)
            except KeyError as e:
                return None
        elif item in self.idmap:
            return self.id(item)
        else:
            self.vcffile.seek(item)
            return Variant.from_str(self.vcffile.readline())

    def __len__(self):
        return len(self.indexmap)

    @property
    def is_phased(self):
        '''
           Samples 100 random variants for phase  
        '''
        variants = [self.rand_variant() for x in range(100)]
        return all([x.is_phased for x in variants])

    @property
    def shape(self):
        return (len(self),len(self.samples))

    @property
    def intervals(self):
        return ["{}:{}".format(chrom,pos) for chrom in self.posmap.keys() for pos in self.posmap[chrom].keys() ]

    def ix(self,index):
        return self[self.indexmap[index]]

    def loc(self,chrom,pos):
        '''
        Access methos to retrieve a variant by its genomic position.
            
        Parameters
        ----------
        chrom : str
            Chromosome name
        pos : int
            Chromosome coordinate
                
        Returns
        -------
            A Variant object at that position
            
        Notes
        -----
        
        '''
        return self.pos(chrom,pos)

    def pos(self,chrom,pos):
        return self[self.posmap[chrom][int(pos)]]

    def id(self,id):
        return self[self.idmap[id]]

    def genos(self,id_list,sample_list,transform=None):
        if not transform:
            transform = lambda x:x
        return pd.DataFrame([self.id(id).genos([self.samples.index(s) for s in sample_list],transform=transform)
            for id in id_list],index=id_list,columns=sample_list)

    def concord_AxiomCalls(self,calls,annot,by_axiomid=False):
        '''
        Compares the genotype calls in the VCF to comparable calls
        in an AxiomCalls object.
            
        Parameters
        ----------
        calls : ponytools.AxiomCalls object
            The object the cross reference here.
                
        Returns
        -------
            A PCC table comparing sample genotypes across objects. 
        '''

        # Grab common samples
        common_samples = set(self.samples).intersection(calls.samples)

        log('Smashing annot and geno together')
        # Smash together calls 
        right_genos = annot._annot.reset_index()[
            ['Probe Set ID','Allele A','Allele B','cust_chr','cust_pos']
        ].set_index('Probe Set ID').join(calls._calls,how='inner').reset_index().set_index(['cust_chr','cust_pos'])

        log('Grabbing vcf variants') 
        # Grab common variants
        vcf_loci = [(chrom,pos) for chrom in self.posmap for pos in self.posmap[chrom]]
        vcf_probeid = [(chrom,pos) for chrom in self.posmap for pos in self.posmap[chrom]]

        common_loci = list(set(vcf_loci).intersection(set(right_genos.index.values)))

        log('Reshaping annot')
        right_genos = right_genos.groupby(level=right_genos.index.names).last()
        right_genos = right_genos.loc[common_loci]
        assert len(common_loci) == len(right_genos)

        # extract series for ref calls
        ref_alleles = right_genos['Allele A']

        log('building vcf -> geno table')
        vcf_vars = [self.pos(chrom,pos).conform(ref_alleles.loc[(chrom,pos)]) for chrom,pos in common_loci]
        assert len(common_loci) == len(vcf_vars)
        vcf_sample_i = [self.sample2index(x) for x in common_samples]
        left_genos = pd.DataFrame(
            [x.alleles(samples_i=vcf_sample_i,transform=Allele.vcf2geno) for x in vcf_vars],
            columns = common_samples,
            index = right_genos.index.values
        )

        results = [('sample','concordance','opp_hom_per','het_per','n')]
        for sample in common_samples:
            # extract nan masks
            nan_mask = np.logical_not(
                np.logical_or(np.isnan(left_genos[sample]), np.isnan(right_genos[sample]))
            )
            n = sum(nan_mask)
            concordance_percent = sum(left_genos[sample][nan_mask] == right_genos[sample][nan_mask])/n
            opposite_hom_percent = sum(abs(left_genos[sample][nan_mask] - right_genos[sample][nan_mask]) == 2)/n
            het_percent = sum(abs(left_genos[sample][nan_mask] - right_genos[sample][nan_mask]) == 1)/n
            results.append((sample,concordance_percent,opposite_hom_percent,het_percent,n))
        return results


    def concord_VCF(self,vcf):
        '''
            Compares a VCF object to another VCF object, comparing calls.        
            
            Parameters
            ----------
            VCF : ponytools.VCF instance
                The VCF instance you want to compare to
                    
            Returns
            -------
            A PCC table comparing sample genotypes across objects.
        '''
        # iterate through the 
        common_samples = set(self.samples).intersection(vcf.samples)
        if len(common_samples) == 0:
            raise ValueError('No common samples between datasets.')
        log('Generating loci lists')
        left_loci = [(chrom,pos) for chrom in self.posmap for pos in self.posmap[chrom]]
        right_loci = [(chrom,pos) for chrom in vcf.posmap for pos in vcf.posmap[chrom]]
        common_loci = set(left_loci).intersection(
            right_loci
        )

        # Right now does not handle multiple variants at a position
        log('Extracting "left-side" variant calls (n={})',len(left_loci))
        left_vars = [self.pos(chrom,pos) for chrom,pos in common_loci]
        left_vars = sorted(left_vars, key=lambda x:vcf.posmap[x.chrom][x.pos])
        log('Extracting "right-side" variant calls (n={})',len(left_loci))
        right_vars = [vcf.pos(v.chrom,v.pos) for v in left_vars]
        # filter out triallelic and non-polymorphic
        log('Filtering out tri-allelic and non-polymorphic SNPs (n={})',len(left_vars))
        left_vars,right_vars = zip(*[ (l,r.conform(l.ref)) for l,r in zip(left_vars,right_vars) if l.is_biallelic and l.is_polymorphic and r.is_biallelic and r.is_polymorphic])

        log('Checking shapes (n={})',len(left_vars))
        # Assert the fuck out of this
        assert len(left_vars) == len(right_vars),\
                "Variant lists not the same length"

        assert all(np.array([x.ref for x in left_vars]) == np.array([x.ref for x in right_vars])),\
                "VCFs are not conformed correctly"
    
        common_loci = [(x.chrom,x.pos) for x in left_vars]

        left_sample_i = [self.sample2index(x) for x in common_samples]
        right_sample_i = [vcf.sample2index(x) for x in common_samples]

        log('Loading genotypes (n={} loci)',len(common_loci))
        # finally load genos
        left_genos = pd.DataFrame(
            [x.alleles(samples_i=left_sample_i,transform=Allele.vcf2geno) for x in left_vars],
            columns = common_samples,
            index = common_loci
        )
        left_genos[left_genos == -1] = np.nan
        right_genos = pd.DataFrame(
            [x.alleles(samples_i=right_sample_i,transform=Allele.vcf2geno) for x in right_vars],
            columns = common_samples,
            index = common_loci
        )
        right_genos[right_genos == -1] = np.nan

        results = [('sample','concordance','opp_hom_per','het_per','n')]
        for sample in common_samples:
            # extract nan masks
            nan_mask = np.logical_not(
                np.logical_or(np.isnan(left_genos[sample]), np.isnan(right_genos[sample]))
            )
            n = sum(nan_mask)
            concordance_percent = sum(left_genos[sample][nan_mask] == right_genos[sample][nan_mask])/n
            opposite_hom_percent = sum(abs(left_genos[sample][nan_mask] - right_genos[sample][nan_mask]) == 2)/n
            het_percent = sum(abs(left_genos[sample][nan_mask] - right_genos[sample][nan_mask]) == 1)/n
            results.append((sample,concordance_percent,opposite_hom_percent,het_percent,n))
        return results
        
    def check_trio(self,offspring,father,mother,return_raw=False):
        ''' this checks the consitency of trio data '''
        raise NotImplementedError()
        consistencies = []
        off_i,fat_i,mot_i = [self.sample2index(x) for x in [offspring,father,mother]]
        for variant in self.iter_variants():
            if not variant.is_biallelic:
                continue
            off,fat,mot = variant.alleles([off_i,fat_i,mot_i])
            consistencies.append(Trio(off,fat,mot).consistent())
        if return_raw:
            return consistencies
        return sum(consistencies)/len(consistencies)

    def to_fam(self,filename,fam_id="VCF"):
        raise NotImplementedError()
        with open(filename,'w') as OUT:
            for sample in self.samples:
                print("{}\t{}\t0\t0\t0\t0".format(fam_id,sample),file=OUT) 
