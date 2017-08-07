import pandas as pd
import numpy as np
import pkg_resources

class MNEc2MAnnot(object):
    def __init__(self): 
        '''
        The AxiomAnnotation as an API for the flat file that comes with
        '''
        # Internal data structures
        self._annot = pd.read_csv(
            pkg_resources.resource_stream('ponytools','data/MNEc2M_Annotation.csv.gz'),
            compression='gzip',
            low_memory=False,
        )
        self._build_BIECIndex()
        # Create the multiIndex
        self._annot.set_index(['MNEcID','AffySNPID','ProbeSetID','BIECID','chrom','pos'],inplace=True)

    def __getitem__(self,item):
        '''
            Try your damndest to fetch a record
        '''
        if isinstance(item,tuple) and len(item) == 2:
            # looks like coordinates
            chrom,pos = item
            return self._annot.xs((str(chrom),int(pos)),level=('chrom','pos'))
        elif item.upper().startswith('MNEC'):
            return self._annot.xs(item,level='MNEcID')
        elif item.startswith('Affx-'):
            return self._annot.xs(item,level='AffySNPID')
        elif item.startswith('AX-'):
            return self._annot.xs(item,level='ProbeSetID')
        elif item.startswith('BIEC-'):
            return self._annot.xs(item,level='BIECID')
        else:
            raise ValueError('Not a valid item')


    def _build_BIECIndex(self):
        '''
            This generates the legacy BIEC IDs from the annotation file
        '''
        self._annot['BIECID'] = [x.split('.')[4] if 'BIEC' in x else np.nan for x in self._annot.MNEcID]


    def in2M(self,chrom,pos):
        '''
        Returns True if coordinates in the MNEc2M Set

        Parameters
        ----------
        chrom: str
            chromosome name
        pos: int
            chromosome position
        '''
        try:
            _ = self[chrom,int(pos)]
            return True
        except KeyError as e:
            return False

    def in670(self,chrom,pos):
        '''
        Returns True if coordinates in the MNEc670k Set

        Parameters
        ----------
        chrom: str
            chromosome name
        pos: int
            chromosome position
        '''
        try:
            _ = self[chrom,int(pos)]
            return all(_.in670 == True)
        except KeyError as e:
            return False

    def getMNEcid(self,item):
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
        x = self[item]
        cust_id = x.index.get_level_values('MNEcID')
        return cust_id[0]

    def getloc(self, item):
        ''' 
            Convert a probe ID to its corresponding genomic chrom and pos.

            Parameters
            ----------
            item: str
                the item you want to look up. Can be the 
                MNEcID, probe id, AffX ID, etc

            Returns
            -------
                a tuple containing chrom and pos

        '''
        x = self[item]
        chrom = x.index.get_level_values('chrom')[0]
        pos   = x.index.get_level_values('pos')[0]
        return (chrom,pos)

    def getprobe(self, item):
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

