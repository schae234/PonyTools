import pandas as pd
from .AxiomAnnot import AxiomAnnot
from .Fasta import Fasta
from .Tools import log,skipno

from .Exceptions import MissingChromosomeError


class AxiomGeno(object):
    def __init__(self):
        self.samples = pd.DataFrame()
        self.CELMap = {}
        self.genos = pd.DataFrame()
        self.annots = AxiomAnnot()
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
        for id,num_seen in [(key,len(list(group))) for key,group in table.groupby(sorted(new_columns))]:
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
                '##source=PonyTools',
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
                variants = self.annots.table[
                    self.annots.table.cust_chr == chrom_name
                ][['cust_chr','cust_pos','Allele A','Allele B']].join(
                    self.genos
                ).sort('cust_pos',ascending=True)
                for id,var in variants.iterrows():
                    chrom,pos,allele_A,allele_B = var[['cust_chr','cust_pos','Allele A','Allele B']] 
                    genotypes = [str(self.genomap[x])+':42' for x in var[self.genos.columns]]
                    if allele_A.upper() not in ['A','C','G','T'] or allele_B.upper() not in ['A','C','G','T']:
                        continue
                    if self.fasta[chrom][pos].upper() == allele_A.upper():
                        ref = allele_A.upper()
                        alt = allele_B.upper()
                        print("\t".join(list(map(str,[chrom,int(pos),id,ref,alt,'.','PASS','SNP=1','GT:GQ']))+genotypes),file=OUT)
                    elif self.fasta[chrom][pos].upper() == allele_B.upper():
                        ref = allele_B.upper()
                        alt = allele_A.upper()
                    else:
                        log("Probe alleles do not match ref for {}({}:{}): ref:{} A:{} B:{}".format(id,chrom,pos,self.fasta[chrom][pos],allele_A,allele_B))
                        continue



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
        for file in geno_files:
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


