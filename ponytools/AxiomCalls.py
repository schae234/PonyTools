import pandas as pd
from collections import defaultdict

from .Tools import skipno,log
from .Variant import Variant
from .Allele import Allele
from .Exceptions import MissingChromosomeError,TriAllelicError
from .VCF import VCF

from ponytools import __version__

class AxiomCalls(object):
    def __init__(self,df):
        self._file_names = []
        self._header = dict()
        self._update_calls(df)


    @staticmethod
    def uniqify(iter):
        'returns a uniqified version of'
        seen = defaultdict(int)
        unique = []
        for item in iter:
            if seen[item] == 0:
                unique.append(item)
            else:
                unique.append('{}_{}'.format(item,seen[item])) 
            seen[item] += 1
        return unique 
        

    @property
    def samples(self):
        # sample names are stored in the columns of calls
        return self._calls.columns

    def _update_calls(self,df):
        '''
        Sets genotype calls in instance from a Pandas df.
            
        Parameters
        ----------
        df : pandas.DataFrame
            A dataframe containing the calls
                
        Returns
        -------
        Self. An updated AxiomCalls intance. 
            
        Notes
        -----
        See AxiomCalls.add_calls for adding genotypes.
        '''
        self._calls = df
        # make sure that sample names are unique
        self._calls.columns = self.uniqify(self._calls.columns)

        

    @classmethod
    def from_file(cls,call_file,sample_file=None):
        calls = pd.read_table(call_file,skiprows=skipno(call_file)).set_index('probeset_id') 
        self = cls(calls)
        self._file_names.append(call_file)
        if sample_file is not None:
            # Replace the column names with sample names
            log('Replacing CELs with Samples names')
            if sample_file.endswith('.xls'):
                samples = pd.read_excel(sample_file)
            elif sample_file.endswith('.xlsx'):
                samples = pd.read_excel(sample_file)
            else:
                samples = pd.read_table(sample_file,sep=',')
            # Remove the spaces from the Column Names
            samples.columns = [x.replace(' ','') for x in samples.columns.values]
            samples = samples[['BestArray','SampleName']].set_index('BestArray').to_dict()
            self._calls.columns = self.uniqify(
                [samples['SampleName'].get(CEL.replace('.CEL',''),CEL) for CEL in self._calls.columns]
            )
        return self
 

    def to_vcf(self,filename,AxiomAnnot,Fasta,conform=True,skip_triallelic=True):
        '''
        Converts the AxiomCall object to a VCF file. Requires an AxiomAnnot
        object so we can assign proper VCF file information.         
            
        Parameters
        ----------
        filename : str
            output filename
        AxiomAnnot : ponytool.AxiomAnnot instance
            Probe information.
        Fasta : ponytools.Fasta instance (default: None)
            Proper reference information. If None, default REF is A Allele.
        conform: bool (default=True)
            Arranges the genotype calls the the REF is the allele found in
            the Fasta object at that locus.
        skip_triallelic : bool (default=True)
            If there is a tri-allelic SNP, just skip it and keep going.
                
        Returns
        -------
            None
            
        '''
        with open(filename,'w') as OUT:
            print("\n".join([
                '##fileformat=VCFv4.0',
                '##source=PonyTools:v{}'.format(__version__),
                '##AxiomCallFiles={}'.format(','.join(self._file_names)),
                '##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">',
                '##INFO=<ID=MNEc,Number=0,Type=Flag,Description="Variant is a MNEc SNP">',
                '##INFO=<ID=MNEcID,Number=1,Type=String,Description="An id the is actually useful">'
            ]),file=OUT)
            # print out chromosomes and their lengths
           
            print("\t".join(['#CHROM','POS','ID','REF','ALT', \
                    'QUAL','FILTER','INFO','FORMAT']+list(self.samples)),file=OUT)
            # reset index and join tables
            variants = AxiomAnnot._annot.reset_index().set_index('Probe Set ID')
            for chrom_name in Fasta.added_order:
                log('Printing {}',chrom_name)
                chrom_vars = variants[variants.cust_chr == chrom_name].join(
                    self._calls,how='inner'
                ).sort_values(by='cust_pos',ascending=True)
                for i,(probeid,var) in enumerate(chrom_vars.iterrows()):
                    if i % 1000 == 0:
                        log("Processed {} variants...".format(i))
                    MNEcID,chrom,pos,allele_A,allele_B = var[
                        ['cust_id','cust_chr','cust_pos','Allele A','Allele B']
                    ] 
                    genotypes = var.loc[self.samples]
                    assert genotypes.shape[0] == len(self.samples)
                    variant = Variant(
                        chrom, pos, probeid, 
                        allele_A, allele_B,
                        '.', 'PASS',info='.',fmt='GT:GQ',
                        genos=["{}:42".format(Allele.geno2vcf(x)) for x in genotypes]
                    )
                    variant.add_info('MNEcID',MNEcID)
                    if variant.is_polymorphic is False:
                        continue
                    if conform is True:
                        try:
                            #Conform to Fasta REF
                            fasta_ref = Fasta[variant.chrom][variant.pos]
                            variant = variant.conform(fasta_ref.upper())
                        except MissingChromosomeError as e:
                            pass
                        except TriAllelicError as e:
                            if skip_triallelic is True:
                                continue
                            else:
                                raise e
                    print(variant,file=OUT)

    


