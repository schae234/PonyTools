import pandas as pd
from .Tools import skipno,log

class AxiomAnnot(object):
    def __init__(self,df,filename=None): 
        '''
        The AxiomAnnotation as an API for the flat file that comes with
        '''
        # Internal data structures
        self._file_names = []
        if filename is not None:
            self._file_names.append(filename) 
        # Build
        self._annot = df.convert_objects()
        self._set_index()

    def ix(self,probeid):
        return self._annot.xs(probeid,level='Probe Set ID')

    def _set_index(self):
        # set up the main internal index
        self._annot.set_index(
            ['cust_chr','cust_pos','cust_id','Probe Set ID'],
            inplace=True
        )

    def probe2MNEcid(self,probeid):
        '''
        Converts between a probe id and the MNEc id.        
            
        Parameters
        ----------
        probeid : str
            The affy probe id
                
        Returns
        -------
        A string containing the MNEc id.
        '''
        x = self._annot.xs(probeid,level='Probe Set ID')
        cust_id = x.index.get_level_values('cust_id')
        return cust_id[0]


    def probe2loc(self, probeid):
        ''' 
            Convert a probe ID to its corresponding genomic chrom and pos.

            Parameters
            ----------
            probeid: str
                the probe you want to look up

            Returns
            -------
                a tuple containing chrom and pos

        '''
        x = self._annot.xs(probeid,level='Probe Set ID')
        chrom = x.index.get_level_values('cust_chr')
        pos   = x.index.get_level_values('cust_pos')
        if len(chrom) != 1:
            raise ValueError(
                "Probe {} has more than one coordinate.".format(probeid)
            )
        return (chrom[0],pos[0])

    def loc2probe(self, chrom, pos, allow_multiple=False):
        '''
            Converts genomic positions to corresponding probe ID

            Parameters
            ----------
                chrom : str
                    chromosome name
                pos : int
                    genetic position
                allow_multiple : bool (default: False)
                    Flag to allow multiple probes to map to locus
                    
            Returns
            -------
                Corresponding probe id
                
            Notes
            -----
        
        '''
        x = self._annot.xs((str(chrom),int(pos)),level=('cust_chr','cust_pos'))
        probe = x.index.get_level_values('Probe Set ID')
        if len(probe) > 1 and False == allow_multiple:
            raise ValueError(
                "Genomic positions have multiple probes, see docstring"
            )
        if len(probe) == 1:
            return probe.values[0]
        else:
            return probe.values


    @classmethod
    def from_file(cls,filename,sep=','):
        # create an empty data structure
        df = pd.read_table(filename,sep=sep,skiprows=skipno(filename))
        self = cls(df, filename=filename)
        return self






