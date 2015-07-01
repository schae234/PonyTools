
from ponytools.Exceptions import TriAllelicError

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
    def chrom(self):
        return self['chrom']

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
        return " ".join(map(str,[self.chrom, self.pos, self.id, self.ref, 
                    self.alt, self.qual, self.filter, self.info, 
                    self.format ]+ self.genotypes[0:5] ))

    def __str__(self):
        '''
            returns a string of variant suitable for VCF printing
        '''
        return "\t".join(map(str,[self.chrom, self.pos, self.id, self.ref, 
                    self.alt, self.qual, self.filter, self.info, 
                    self.format ]+ self.genotypes ))

    def alt_freq(self,samples_i=None,max_missing=0.3):
        '''
            Computes the ALTERNATE allele frequency.

            Parameters
            ----------
            samples_i : str iterable (individual ids)
                If not None, will compute MAF for subset of individuals.

        '''
        genos = self.genos()
        if samples_i is not None:
            genos = [genos[x] for x in samples_i]
        num_mis = sum([x.count('.') for x in genos])
        if (num_mis / (2*len(genos))) > max_missing:
            return None
        num_ref = sum([x.count('0') for x in genos])
        num_alt = sum([x.count('1') for x in genos])
        n = (2*len(genos)) - num_mis
        return num_alt/n

    def heterozygosity(self,samples_i=None,max_missing=0.3):
        '''
            Compute the heterozygositey of the SNP
            for indivduials indices in samples_i

            parameters
            ----------
            samples_i : int iterable (individual ids)
                If not None, will computer heterozygosity for 
                subset of individuals
        '''
        alt_freq = self.alt_freq(samples_i=samples_i,max_missing=max_missing)
        if alt_freq is None:
            return None
        ref_freq = 1 - alt_freq
        return 2 * alt_freq * ref_freq


    def genos(self,samples_i=None,samples_id=None,format_field='GT',transform=None):
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



    def conform(self,reference_genotype,GT_index=0):
        '''
            Conforms a variant to a reference genotype.
            This happens IN PLACE.

            Parameters
            ----------
            reference_genotype: str
                must be current alt or ref genotype. This
                becomes the new reference.

        '''
        if self.alt == self.ref and self.ref != reference_genotype:
            # Sometimes PLINK will assign both the alt and reference
            # to the same thing, just make everyone alternate
            self.ref = reference_genotype
            # Make the translation table
            trans = str.maketrans(
                {'0':'1',
                 '1':'1'}
            )
            def conform_genotype(genotype):
                # split out the GT 
                fields = genotype.split(':')
                fields[GT_index] = str.translate(fields[GT_index],trans)
                return ':'.join(fields)
            self.genotypes = [conform_genotype(g) for g in self.genotypes]

        elif reference_genotype not in (self.alt,self.ref):
            raise TriAllelicError((
                'Cannot conform to genotype that is not one of'
                'the reference or alt {} != ({},{})'.format(
                    reference_genotype,
                    self.ref,
                    self.alt
                )
            ))
        elif not self.biallelic:
            raise TriAllelicError('cannot conform triallelic SNP: {}'.format(self.id))
        elif self.ref == reference_genotype:
            return
        else:
            # Conform
            self.alt,self.ref = self.ref,self.alt
            # Make the translation table
            trans = str.maketrans(
                {'0':'1',
                 '1':'0'}
            )
            def conform_genotype(genotype):
                # split out the GT 
                fields = genotype.split(':')
                fields[GT_index] = str.translate(fields[GT_index],trans)
                return ':'.join(fields)
            self.genotypes = [conform_genotype(g) for g in self.genotypes]
            
            
             
        
 
