class Variant(object):
    genomap = { # shared genomap
        './.' : -1, '0/0' : 0, '0/1' : 1, '1/0' : 1, '1/1' : 2,
        -1 : -1, 0 : 0, 1 : 1, 2 : 2,
        '.|.' : -1, '0|0' : 0, '0|1' : 1, '1|0' : 1, '1|1' : 2,
    }
    def __init__(self,string):
        # make this lightweight since we are going to be using a  these
        fields = string.strip().split()
        self.__dict__.update(zip(['chrom','pos','id','ref','alt','qual','filter','info','format'],fields[0:9]))
        self.genotypes = fields[9:]

    @property
    def phased(self):
        ''' return bool on phase state '''
        if '|' in self.genotypes[0]:
            return True
        else:
            return False
    @property
    def biallelic(self):
        ''' return bool true if variant is biallelic '''
        if ',' in self.alt:
            return False
        else:
            return True
    @property
    def pos(self):
        return int(self.__dict__['pos'])
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
        return "\n".join(map(str,[self.chrom, self.pos, self.id, self.ref, 
                    self.alt, self.qual, self.filter, self.info, 
                    self.format ]+ self.genotypes))
    def genos(self,samples_i=None,format_field='GT',transform=None):
        ''' return alleles *in order* for samples '''
        if not transform:
            transform = lambda x:x
        gt_index = self.format.split(":").index(format_field)
        if samples_i:
            return [ transform(self.genotypes[sample_i].split(':')[gt_index]) for sample_i in samples_i]
        else:
            return [transform(geno.split(':')[gt_index]) for geno in self.genotypes]
    def r2(self,variant,samples_i,samples_j):
        ''' returns the r2 value with another variant '''
        # find common individuals
        geno1 = np.array(self.genos(transform=Allele.vcf2geno,samples_i=samples_i))
        geno2 = np.array(variant.genos(transform=Allele.vcf2geno,samples_i=samples_j))
        non_missing = (geno1!=-1)&(geno2!=-1)
        return (scipy.stats.pearsonr(geno1[non_missing],geno2[non_missing]))[0]**2
    
    def __sub__(self,variant):
        ''' returns the distance between SNPs '''
        if self.chrom != variant.chrom:
            return np.inf
        else:
            return abs(int(self.pos) - int(variant.pos))
        
 
