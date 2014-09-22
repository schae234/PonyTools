#!/usr/bin/env python3
import sys
import os
import scipy.stats
from collections import defaultdict
from itertools import chain,groupby
import pandas as pd
import gzip
import pickle
import numpy as np

def skipno(filename,startswith="#"):
    with open(filename,'r') as IN:
        skip = 0
        while IN.readline().startswith(startswith):
            skip += 1
        return skip

def log(message,*formatting):
    print(message.format(*formatting),file=sys.stderr)


def SpecialSNPs():
    special_snps = pd.DataFrame()
    special_path = '/project/mccuelab/rob/MNEc2M/special_snps.txt'
    for i,design in enumerate(['A','B','C']):
        annot = "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Annotation.v1/CsvAnnotationFile_{}.v1.txt".format(design,design,i+1)
        special_snps = pd.concat([special_snps,pd.read_csv(annot,skiprows=skipno(annot),sep=',',quotechar='"')])
    special_snps = special_snps[[x in ['chrX','chrY'] for x in special_snps.cust_chr]]
    special_snps = special_snps[["Affy SNP ID","cust_chr"]]
    special_snps.columns = ['probeset_id','chr']
    special_snps['copy_female'] = [2 if chr == 'chrX' else 0  for chr in special_snps.chr ]
    special_snps['copy_male'] = 1
    special_snps.to_csv(special_path,sep="\t",index=False)
 
def read_vcf(filename):
    mapping = []
    with open(filename,'r') as IN:
        for line in IN:
            if line.startswith('#'):
                continue
            chrom,pos,id,ref,alt,qual,filter,info,format,*genos = line.strip().split()
            info = info.replace(';','=').split('=')
            info = dict(zip(info[::2],info[1::2]))
            info['chrom'] = chrom
            info['pos'] = pos
            info['qual'] = qual
            mapping.append(info)
    return mapping

def read_genotypes(vcf,mapfile,outfile="WGSCalls.tsv"):
    import re
    outfile = open(outfile,'w')
    geno_pattern = re.compile(":.*")
    geno_map = {
        './.' : '-1',
        '0/0' : '0',
        '0/1' : '1',
        '1/0' : '1',
        '1/1' : '2'
    }
    snpmap = dict()
    with open(map_file,'r') as MAP:
        for i,line in enumerate(MAP):
            chrom,pos,cm,id = line.strip().split('\t')
            snpmap[id] = i 
    with open(vcf,'r') as IN:
        for i,line in enumerate(IN):
            if i % 100000 == 0:
                print("On line {}".format(i))
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                chrom,pos,id,ref,alt,qual,filter,info,format,*samples = line.strip().split()
            else:
                chrom,pos,id,ref,alt,qual,filter,info,format,*genos = line.strip().split()
                # Check if we have something worth printing out!
                if (chrom,int(pos)) in idmap:
                    # Cut out the genotpyes
                    probeid = idmap[(chrom,int(pos))]
                    for ind,geno in zip(samples,genos):
                        geno = re.sub(geno_pattern,'',geno)
                        print("{}\t{}\t{}".format(probeid,ind,geno_map[geno]),file=outfile)
                    

def long_to_ped(long_file_list,map_file,outped="Output.ped"):
    # read map file so we know dims of eventual matrix
    log("Reading in MAp file")
    snpmap = dict()
    with open(map_file,'r') as MAP:
        for i,line in enumerate(MAP):
            chrom,pos,cm,id = line.strip().split('\t')
            snpmap[id] = i
    log("Read in {} SNPs",len(snpmap))
    indmap = defaultdict(lambda : [-1]*len(snpmap))
    # Deal with the long files
    for long_file in long_file_list:
        log("Reading in long file: {}",long_file)
        with open(long_file,'r') as IN:
            for i,line in enumerate(IN):
                if i % 10000000 == 0:
                    log("On Line {}",i)
                snpid,indid,call = line.strip().split('\t')
                # Skip the snps not in our map file
                if snpid not in snpmap:
                    continue
                else:
                    indmap[indid][snpmap[snpid]] = int(call)
    # print out ped file
    genomap = {-1 : '0 0', 0 : '1 1', 1 : '1 2', 2 : '2 2'}
    with open(outped,'w') as OUTPED:
        for ind,genos in indmap.items():
            log("Writing genotypes for {}",ind)
            print("0\t{}\t0\t0\t0\t0\t{}".format(ind,"\t".join([genomap[x] for x in genos])),file=OUTPED)
         
 
def HW_chi(n_AA,n_AB,n_BB):
    n = sum([n_AA,n_AB,n_BB])
    p = (2*n_AA+n_AB)/(2*n)
    q = 1-p
    exp_AA = (p**2)*n
    exp_AB = 2*p*q*n
    exp_BB = (q**2)*n
    return scipy.stats.chisquare([n_AA,n_AB,n_BB],[exp_AA,exp_AB,exp_BB],ddof=1)

def VCFFilterByMap(vcf_file,map_file,output_file):
    include = set()
    num_found = 0
    log('reading in map file: {}',map_file)
    with open(map_file,'r') as MAP:
        for line in MAP:
            chrom,id,cm,pos = line.strip().split()
            if not chrom.startswith('chr'):
                chrom = 'chr'+chrom
            include.add((chrom,pos))
    log('reading in vcf file: {}',vcf_file)
    with open(vcf_file,'r') as VCF:
        with open(output_file,'w') as OUT:
            for i,line in enumerate(VCF):
                if i % 100000 == 0:
                    log("on line {}, found {}",i,num_found)
                line = line.strip()
                if line.startswith('#'):
                    print(line,file=OUT)
                else:
                    chrom,pos,*junk = line.split()
                    if (chrom,pos) in include:
                        num_found += 1
                        print(line,file=OUT)
            

def geno2VCF(geno_file,annotation_file,sample_file,reference_file,outfile=None):
    """ This function converts the affymetrix genotype file to a standard VCF file. 
        It requires that you also provide the appropriate annotattion file bundled
        with the R1 package.
    """
    # Read in the Fasta File
    refgen = dict()
    with open(reference_file,'r') as FASTA:
        log("Reading in FASTA: {}",reference_file)
        cur_chrom = None
        for line in FASTA:
            line = line.strip()
            if line.startswith('>'):
                cur_chrom = line.replace('>','').split()[0]
                log("Found a new chromosome: {}",cur_chrom)
                refgen[cur_chrom] = list([None]) # positions start at 1 in genomes
            else:
                refgen[cur_chrom].extend(line)
    log("Read in {} chromosomes: [{}]",len(refgen),",".join(refgen.keys()))
    # Create vcf filename from geno_file
    if not outfile:
        outfile = geno_file+'.vcf'
    # Read in Samples to create CEL map:
    log("Reading in sample file: {}",sample_file)
    samples = pd.read_table(sample_file,sep=',')
    CELMap = {}
    for CEL,sample in samples[["BestArray",'SampleName']].itertuples(index=False):
        CELMap[CEL] = sample
    uniq_samples = set(CELMap.values())
    log("Found {} samples: [{}]",len(uniq_samples),",".join(uniq_samples))
    # Read in annotations for 
    log("Reading in annotations: {}",annotation_file)
    annotations = pd.read_table(annotation_file,skiprows=skipno(annotation_file),sep=',')
    # We want to be able to extract probes by their ids
    annotations.index = annotations['Probe Set ID']
    genomap = {'-1' : './.', '0' : '0/0', '1': '0/1', '2' : '1/1'}
    with open(outfile,'w') as VCF:
        print("\n".join([
            '##fileformat=VCFv4.0',
            '##source=PonyTools'
            '##geno_file={}'.format(geno_file),
            '##annotation_file={}'.format(annotation_file),
            '##sample_file={}'.format(sample_file),
            '##reference_file={}'.format(reference_file),
            '##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">',
            '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">',
            '##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">'
        ]),file=VCF)
        with open(geno_file,'r') as GENO:
            for i,line in enumerate(GENO):
                if i % 100000 == 0:
                    log("processed {} probes",i)
                if line.startswith('#'):
                    continue
                elif line.startswith('probeset_id'):
                    samples = [CELMap[cel] for cel in line.strip().split()[1:]]
                    print("\t".join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+samples),file=VCF)
                else:
                    probe,*genotypes = line.strip().split()
                    genotypes = [str(genomap[x])+':42' for x in genotypes]
                    chrom,pos,id,allele_A,allele_B = annotations.ix[probe][['cust_chr','cust_pos','Probe Set ID','Allele A','Allele B']]
                    # figure out which allele is the reference
                    if chrom not in refgen:
                        continue
                    # skip indels:
                    if allele_A.upper() not in ['A','C','G','T'] or allele_B.upper() not in ['A','C','G','T']:
                        continue
                    if refgen[chrom][pos].upper() == allele_A.upper():
                        ref = allele_A.upper()
                        alt = allele_B.upper()
                    elif refgen[chrom][pos].upper() == allele_B.upper():
                        ref = allele_B.upper()
                        alt = allele_A.upper()
                    else:
                        log("Probe alleles do not match ref for {}({}:{}): ref:{} A:{} B:{}".format(probe,chrom,pos,refgen[chrom][pos],allele_A,allele_B))
                        continue
                    print("\t".join(list(map(str,[chrom,pos,id,ref,alt,'.','PASS','SNP=1','GT:GQ']))+genotypes),file=VCF)

def sortVCF(vcf_file,fasta_file,temp_dir="/tmp",out="sorted.vcf"):
    headers = list()
    variants = list()
    cur_byte = 0
    chroms = list() 
    temps  = dict()
    log("Sorting {}",vcf_file)
    # Get the chromosome order
    with open(fasta_file,'r') as FASTA:
        for line in FASTA:
            if line.startswith('>'):
                chrom,*info = line.strip().lstrip('>').split()
                log("Found chromosome {}",chrom)
                chroms.append(chrom)
                temps[chrom] = open(os.path.join(temp_dir,chrom+'.tmp'),'w')
    # Get headers and extract positions with file byte offsets
    log("Reading in VCF: {}",vcf_file)
    with open(vcf_file,'r') as VCF:
        for i,line in enumerate(VCF):
            line = line.strip()
            if line.startswith("#"):
                headers.append(line)
            else:
                chrom,pos,*junk = line.split()
                print(line,file=temps[chrom])
    # close all temp files
    for key,val in temps.items():
        log("Closing tmp file: {}",key)
        val.close()
    log("soring chroms")
    with open(out,'w') as OUT:
        # print headers
        print("\n".join(headers),file=OUT)
        for chrom in chroms:
            # read in that chroms bullshit
            with open(os.path.join(temp_dir,chrom+'.tmp'),'r') as CHROM:
                variants =  CHROM.readlines()
                # sort by position
                variants.sort(key=lambda x: int(x.split()[1]))
                log("printing chrom {}",chrom)
                print("".join(variants),file=OUT,end="")
                os.remove(os.path.join(temp_dir,chrom+'.tmp'))
            
           
def phasedVCF2FastTagger(phased_VCF,fasta,annot,outprefix=None):
    if not outprefix:
        # use the base name of phased_VCF
        outprefix = os.path.basename(phased_VCF).replace('vcf','fasttagger')
    if 'gz' in phased_VCF:
        vcf_open = gzip.open
        mode = 'rb'
    else:
        vcf_open = open   
        mode = 'r'
    with vcf_open(phased_VCF,mode) as VCF, \
         open(outprefix+'.data_maf','w') as MAP, \
         open(outprefix+'.data_matrix','w') as MAT:
        for line in VCF:
            line = line.decode()
            if line.startswith('#'):
                continue
            chrom,pos,id,ref,alt,qual,filters,info,formats,*genos = line.strip().split()
            info = info.replace('=',';').split(';')
            info = dict(zip(info[0::2],info[1::2]))
            GT_index = formats.split(":").index("GT")
            if ref != fasta[chrom][pos]:
                log('reference bases do not match: {}'," ".join(line.split()[0:9]))
            elif ',' in alt:
                if ',' in id:
                    # Find out if is VIP SNP, skip dual probe SNPs
                    continue
                mafs = [float(x) for x in info['AF'].split(',')]
                # MAFs below 0.005 are probably sequencing errors
                if min(mafs) < 0.005:
                    alt = alt.split(',')[mafs.index(max(mafs))]
                    info['AF'] = max(mafs)
                    genos = [geno.replace('2','.') for geno in genos] # replace the third allele with missing
                else:
                    log("{} is tri-allelic!: {}",id," ".join([chrom,pos,id,ref,alt,qual,filters,";".join(["{}={}".format(k,v) for k,v in info.items()]),formats]))
                    for probe in id.split(','):
                        log("\t{}"," ".join(map(str,annot.table.ix[probe])))
                        log("\t{}",' '.join([geno.split(':')[GT_index] for geno in genos]))
                    continue
            # fastTagger is peculiar about its map files
            if float(info['AF']) > 0.5:
                minor_allele = ref
                major_allele = alt
                maf = 1.0 - float(info['AF'])
            else:
                minor_allele = alt
                major_allele = ref
                maf = float(info['AF'])
            # print line to MAP file
            print("{}\t{}\t{}\t{}\t{}".format(id,pos,major_allele,minor_allele,maf),file=MAP)
            # print corresponding line to matrix file
            print("\t".join(chain(*[geno.split(":")[GT_index].split('|') for geno in genos])),file=MAT)
           


class Variant(object):
    def __init__(self,string):
        # make this lightweight since we are going to be using a  these
        fields = string.strip().split()
        self.__dict__.update(zip(['chrom','pos','id','ref','alt','qual','filter','info','format'],fields[0:9]))
        self.genotypes = fields[9:]
    def __getitem__(self,item):
        try:
            return self.__dict__[item]
        except KeyError as e:
            pass
        try:
            return self.genos[item]
        except Exception as e:
            pass
        raise Exception("Value not found")
    def __repr__(self):
        return "\n".join([self.chrom, self.pos, self.id, self.ref, 
                    self.alt, self.qual, self.filter, self.info, 
                    self.format ]+ self.genotypes)
    def genos(self,samples_i,format_field='GT',transform=None):
        ''' return alleles *in order* for samples '''
        if not transform:
            transform = lambda x:x
        index = self.format.split(":").index(format_field)
        return [ transform(self.genotypes[sample_i].split(':')[index]) for sample_i in samples_i]
      
class Trio(object):
    genomap = { # shared genomap
        './.' : -1,
        '0/0' : 0,
        '0/1' : 1,
        '1/0' : 1,
        '1/1' : 2,
        -1 : -1,
        0 : 0,
        1 : 1,
        2 : 2     
    }
    def __init__(self,offspring=None,father=None,mother=None):
        self.off = self.genomap[offspring]
        self.fat = self.genomap[father]
        self.mot = self.genomap[mother]

    def consistent(self):
        if self.off == -1: # give benefit of doubt
            consistent = True
        elif self.off == 0 and (self.mot == 0 or self.mot ==1) and (self.fat == 0 or self.fat == 1):
            # both self.father and self.mother must have 0 allele
            consistent = True
        elif self.off == 1 and (not (self.mot == 0 and self.fat == 0) or not (self.mot == 2 and self.fat == 2)):
            # both self.father and self.mother cannot be homozygous
            consistent = True
        elif self.off == 2 and (self.mot == 1 or self.mot ==2) and (self.fat == 1 or self.fat == 2): 
            # both self.father and self.mother must have 1 allele
            consistent = True
        elif self.off == 0 and ((self.mot == -1 and self.fat != 2) or (self.mot != 2 and self.fat == -1)):
            # if one parent is missing, the other must not be opposite homozygous
            consistent = True
        elif self.off == 2 and ((self.mot == -1 and self.fat != 0) or (self.mot != 0 and self.fat == -1)):
            # if one parent is missing, the other must not be opposite homozygous
            consistent = True
        elif self.off == 1 and (self.mot == -1 or self.fat == -1):
            # when self.offspring is homozygous, all bets are self.off
            consistent = True
        else:
            consistent = False
        return consistent

class Allele(object):
    def __init__(self):
        pass    
    @staticmethod
    def vcf2geno(allele_code):
        return {
            './.' : -1,
            '0/0' : 0,
            '0/1' : 1,
            '1/0' : 1,
            '1/1' : 2,
        }[allele_code]
    @staticmethod
    def geno2vcf(allele_code):
        return {
            -1 :'./.',
            0  : '0/0',
            1  :'0/1',
            1  :'1/0',
            2  :'1/1',
        }[allele_code]

        
            
class VCF(object):
    def __init__(self,vcffile):
        self.vcffile = vcffile
        self.idmap = {}
        self.posmap = defaultdict(dict)
        self.samples = []
        self.index()
    def index(self):
        if not os.path.exists(self.vcffile+'.pdx'):
            # create the index file
            log('Index file does not exist, indexing now')
            cur_byte = 0
            with open(self.vcffile,'r') as IN:
                for line in IN:
                    if line.startswith('##'):
                        pass
                    elif line.startswith("#"):
                        self.samples = line.strip().split()[9:]
                    else:
                        chrom,pos,ids,*junk = line.strip().split()
                        pos = int(pos)
                        for id in ids.split(','): #sometimes there are multiple ids...
                            self.idmap[id] = cur_byte
                        self.posmap[chrom][pos] = cur_byte
                    cur_byte += len(line)
            # pickle index file for later use
            pickle.dump((self.idmap,self.posmap,self.samples),open(self.vcffile+".pdx",'wb'))
        else:
            # read the index file
            self.idmap,self.posmap,self.samples = pickle.load(open(self.vcffile+'.pdx','rb')) 

    def iter_variants(self):
        ''' returns variant generator, for iteration. Should be memory efficient '''
        return (Variant(line) for line in open(self.vcffile,'r') if not line.startswith('#'))

    def sample2index(self,sample):
        ''' return the genotype index for sample '''
        return self.samples.index(sample)

    def __getitem__(self,byte):
        with open(self.vcffile,'r') as VCF:
            VCF.seek(byte)
            return Variant(VCF.readline())

    def pos(self,chrom,pos):
        return self[self.posmap[chr][pos]]

    def id(self,id):
        return self[self.idmap[id]]

    def genos(self,id_list,sample_list,transform=None):
        if not transform:
            transform = lambda x:x
        return pd.DataFrame([self.id(id).genos([self.samples.index(s) for s in sample_list],transform=transform)
            for id in id_list],index=id_list,columns=sample_list)

    def check_trio(self,offspring,father,mother,return_raw=False):
        ''' this checks the consitency of trio data '''
        consistencies = []
        off_i,fat_i,mot_i = [self.sample2index(x) for x in [offspring,father,mother]]
        for variant in self.iter_variants():
            off,fat,mot = variant.genos([off_i,fat_i,mot_i])
            consistencies.append(Trio(off,fat,mot).consistent())
        if return_raw:
            return consistencies
        return sum(consistencies)/len(consistencies)

    
class AxiomGenos(object):
    def __init__(self):
        self.samples = pd.DataFrame()
        self.CELMap = {}
        self.genos = pd.DataFrame()
        self.annots = annotation()
        self.fasta = Fasta()
        self.genomap = {
            -1 : './.',
            0  : '0/0',
            1  : '0/1',
            2  : '1/1'
        }
    def add_sample_table(self,sample_table):
        log("Adding samples")
        # Get rid of spaces in column headers
        sample_table.columns = [x.replace(' ','') for x in sample_table.columns]
        sample_table.index = sample_table['SampleName']
        for CEL,sample in sample_table[['BestArray','SampleName']].itertuples(index=False):
            self.CELMap[CEL] = sample
        self.samples = pd.concat([self.samples,sample_table])
    def CEL2id(self,CEL):
        ''' returns the id for a CEL, if available, otherwise CEL'''
        if CEL in self.CELMap:
            return self.CELMap[CEL]
        elif CEL.replace('.CEL','') in self.CELMap:
            return self.CELMap[CEL.replace('.CEL','')]
        else:
            return CEL
    def geno2AB(self,geno):
        return self.genomap[int(geno)]
    def add_design(self,table):
        ''' This method is for adding additional probes (i.e. a seperate design) 
            to an existing set of samples. It 'appends' probes to the end of the 
            matrix'''
        table.columns  = [self.CEL2id(CEL) for CEL in table.columns] 
        new_columns = list(table.columns)
        # take care of duplicate columns:
        for id,num_seen in [(key,len(list(group))) for key,group in groupby(sorted(new_columns))]:
            if num_seen > 1:
                for i in range(1,num_seen+1):
                    new_columns[new_columns.index(id)] = "{}_{}".format(id,i) 
        log("Colnames set from :")
        log(" [{} ... ]"," ".join(list(table.columns)[0:10]))
        table.columns = new_columns 
        log(" [{} ... ]"," ".join(list(table.columns)[0:10]))
        self.genos = self.genos.append(table)
        self.genos.fillna(-1)

    def add_samples(self,table):
        ''' This method is for adding a table who has the SAME probes,
            but for different samples. i.e. two different experiments. '''
        # join the samples
        self.genos = self.genos.join(table,how='outer')
        self.genos.fillna(-1)
    def add_annotation(self,annot):
        ''' integrates an annotations object '''
        self.annots = self.annots + annot

    def check_trio(self,offspring,father,mother,return_raw=False):
        ''' this checks the consitency of trio data '''
        consistencies = []
        for probe,off,fat,mot in self.genos[[offspring,father,mother]].itertuples():
            off,fat,mot = [int(x) for x in [off,fat,mot]]
            if off == -1: # give benefit of doubt
                consistent = True
            elif off == 0 and (mot == 0 or mot ==1) and (fat == 0 or fat == 1):
                # both father and mother must have 0 allele
                consistent = True
            elif off == 1 and (not (mot == 0 and fat == 0) or not (mot == 2 and fat == 2)):
                # both father and mother cannot be homozygous
                consistent = True
            elif off == 2 and (mot == 1 or mot ==2) and (fat == 1 or fat == 2): 
                # both father and mother must have 1 allele
                consistent = True
            elif off == 0 and ((mot == -1 and fat != 2) or (mot != 2 and fat == -1)):
                # if one parent is missing, the other must not be opposite homozygous
                consistent = True
            elif off == 2 and ((mot == -1 and fat != 0) or (mot != 0 and fat == -1)):
                # if one parent is missing, the other must not be opposite homozygous
                consistent = True
            elif off == 1 and (mot == -1 or fat == -1):
                # when offspring is homozygous, all bets are off
                consistent = True
            else:
                consistent = False
            consistencies.append(consistent)
        if return_raw:
            return consistencies
        return sum(consistencies)/len(consistencies)

    def to_vcf(self,filename):
        # require a fasta, annotations and genotypes
        with open(filename,'w') as OUT:
            print("\n".join([
                '##fileformat=VCFv4.0',
                '##source=PonyTools'
                '##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">',
                '##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">'
            ]),file=OUT)
            print("\t".join(['#CHROM','POS','ID','REF','ALT', \
                    'QUAL','FILTER','INFO','FORMAT']+list(self.genos.columns)),file=OUT)
            # Print the variants in sorted order by fasta chromosome position then by position 
            # annotate the probe table
            log('joining annotations')
            for chrom_name in self.fasta.added_order:
                log('Printing {}',chrom_name)
                # Retrieve the probes associated with that chromosome
                self.annots.table[self.annots.table.cust_chr == chrom_name]
                variants = self.annots.table[self.annots.table.cust_chr == chrom_name][['cust_chr','cust_pos','Allele A','Allele B']].join(self.genos).sort('cust_pos',ascending=True)
                for id,var in variants.iterrows():
                    chrom,pos,allele_A,allele_B = var[['cust_chr','cust_pos','Allele A','Allele B']] 
                    genotypes = [str(self.genomap[x])+':42' for x in var[self.genos.columns]]
                    if allele_A.upper() not in ['A','C','G','T'] or allele_B.upper() not in ['A','C','G','T']:
                        continue
                    if self.fasta[chrom][pos].upper() == allele_A.upper():
                        ref = allele_A.upper()
                        alt = allele_B.upper()
                    elif self.fasta[chrom][pos].upper() == allele_B.upper():
                        ref = allele_B.upper()
                        alt = allele_A.upper()
                    else:
                        log("Probe alleles do not match ref for {}({}:{}): ref:{} A:{} B:{}".format(id,chrom,pos,self.fasta[chrom][pos],allele_A,allele_B))
                        continue
                    print("\t".join(list(map(str,[chrom,int(pos),id,ref,alt,'.','PASS','SNP=1','GT:GQ']))+genotypes),file=OUT)

    def get_pos(self,chr,pos):
        id = self.annots.table[(self.annots.table.cust_chr == str(chr)) & (self.annots.table.cust_pos == int(pos))]["Probe Set ID"][0]
        return self.genos.ix[id]

    def __getitem__(self,item): 
        try:
            ids,probes = item
            return self.genos[ids].ix[probes]
        except Exception as e:
            pass
        
            
    @classmethod
    def from_files(cls,sample_file=None,geno_files=None,fasta_file=None,annot_files=None):
        cls = cls()
        cls.add_sample_table(pd.read_table(sample_file,sep=','))
        for file in files:
            log('adding table: {}',file)
            table = pd.read_table(file,skiprows=skipno(file))
            table.index = table['probeset_id']
            # Check for duplicate columns
            del table['probeset_id']
            cls.add_design(table)
        return cls
    @classmethod
    def from_copy(cls,self):
        cls = cls()
        cls.__dict__ = self.__dict__

class annotation(object):
    def __init__(self):
        self.table = pd.DataFrame()
    def add_table(self,table):
        log("adding new table")
        self.table = pd.concat([self.table,table])
    def __getitem__(self,item):
        return self.table.ix[item]
    def __add__(self,annot):
        return self.add_table(annot.table)

    @classmethod
    def from_files(cls,*files):
        cls = cls()
        for file in files:
            table = pd.read_table(file,skiprows=skipno(file),sep=',')
            table.index = table['Probe Set ID']
            cls.add_table(table)
        return cls

class chromosome(object):
    def __init__(self,seq):
        self.seq = str(seq)
    def __getitem__(self,pos):
        #chromosomes start at 1, python strings start at 0
        return self.seq[int(pos)-1]
       
class Fasta(object):
    def __init__(self):
        self.added_order = []
        self.chroms = {}
        self.attributes = defaultdict(list)

    def add_chrom(self,chrom_name,sequence):
        log("Adding new chromosome: {}",chrom_name)
        self.added_order.append(chrom_name)
        self.chroms[chrom_name] = sequence
    def add_attribute(self,chrom_name,attr):
        self.attributes[chrom_name].append(attr)

    def __getitem__(self,chrom_name):
        return self.chroms[chrom_name]

    @classmethod
    def from_file(cls,fasta_file):
        fasta = cls()
        with open(fasta_file,'r') as IN: 
            cur_chrom = None
            cur_seqs = []
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    if cur_chrom:
                        fasta.add_chrom(cur_chrom,chromosome("".join(cur_seqs)))
                    name,*attrs = line.lstrip('>').split()
                    cur_chrom = name
                    cur_seqs = []
                    for attr in attrs:
                        fasta.add_attribute(attr)
                else:
                    cur_seqs.append(line)
            # Add the last chromosome
            fasta.add_chrom(cur_chrom,chromosome("".join(cur_seqs)))
        return fasta
