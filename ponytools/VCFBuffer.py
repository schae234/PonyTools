#!/usr/bin/python

import sys
import gzip

class VCF(object):
    def __init__(self,filename):
        pass

class Genotype(object):
    def __init__(self,format,fields):
        for key,val in zip(format,fields):
            setattr(self,key,val)
    def get_alleles(self,assume_ref=True):
        missing_value = "-1"
        if assume_ref:
            missing_value = "0"
        alleles = self.GT.split('/')
        return [int(x) for x in self.GT.replace('.',missing_value).split('/')]
    def __str__(self):
        return ":".join(self.__dict__.itervalues())


class Variant(VCF):
    def __init__(self,fields):
        self.chrom = fields[0]
        self.pos   = int(fields[1])
        self.id    = fields[2]
        self.ref   = fields[3]
        self.alt   = fields[4]
        self.qual  = self.process_qual(fields[5])
        self.filters = fields[6].split(';')
        self.info  = self.process_info(fields[7].split(';'))
        self.format = fields[8].split(':')
        #self.genotypes = self.process_genotypes(fields[9:])
        self.genotypes = fields[9:]

    @classmethod
    def from_str(cls,string,sep="/t"):
        string = string.strip()
        return cls(string.split(sep=sep))
        

    def add_filter(self,filter):
        self.filters.append(str(filter))

    def process_qual(self,field):
        if field == '.':
            return '.'
        else:
            return float(field)

    def process_info(self,info_list):
        info_dict = {}
        for x in info_list:
            if '=' in x:
                key,val = x.split("=")
                info_dict[key] = val
            else:
                info_dict[x] = '1'
        info_dict['QUAL'] = self.qual
        return info_dict 

    def process_genotypes(self,list_of_genotypes):
        genotypes = []
        for el in list_of_genotypes:
            genotypes.append(Genotype(self.format,el.split(":")))
        return genotypes

    @property
    def MNEcId(self,tag='PC',gene_build='2'):
        return "MNEc.{}.{}.{}.{}".format(
            gene_build,
            self.chrom.replace('chr',''),
            self.pos,
            tag
        ) 

    @property
    def chr_type(self):
        if 'mito' in self.chrom.lower():
            return 'mitochondrial'
        elif 'un' in self.chrom.lower():
            return 'pseudo-autosomal'
        else:
            try:
                float(self.chrom.replace('chr',''))
                return 'autosomal'
            except ValueError:
                return 'sex'
            

    def __getitem__(self,key):
        return self.fields[key]

    def __str__(self):
        return "\t".join([
            self.chrom,
            str(self.pos), 
            self.id, 
            self.ref, 
            self.alt, 
            str(self.qual), 
            ";".join(self.filters),
            ';'.join(['='.join([str(k),str(v)]) for k,v in self.info.iteritems()]),
            ":".join(self.format),
            "\t".join([ str(g) for g in self.genotypes  ])])

    def __repr__(self):
        return '[{},{},{}...]'.format(self.chrom,self.pos,self.id)

    def get_alleles(self):
        lol = [x.get_alleles() for x in self.genotypes]
        return [item for sublist in lol for item in sublist]


class VCFBuffer(object):
    ''' 
        VCF class which allows you to compare variant entries to
        other entries within a buffer.

        [ [], L1  ]

        '''

    def __init__(self,filename,buffer_function = 'max_size',buffer_size = 10):
        # Set up the buffer clearing function
        self.is_buffer_full = getattr(self,buffer_function)
        # Set up class variable
        self.max_buffer_size = buffer_size # in variants
        self.window_size = buffer_size # in bp
        # Set up class data structures
        self.header_buffer = []
        self.buffer = []
        self.current = 0
        # Set up some class files
        self.log = sys.stderr 
        self.file = self.open(filename,'r')
        # Set up the initial buffer
        self.staged = self.nextvar()
        self.rebuffer()

    @property
    def previ(self):
        if self.current == 0:
            return 0
        else:
            return self.current-1

    @property
    def nexti(self):
        if self.buffer[self.current+1] == []:
            return self.current
        else:
            return self.current + 1

    def nextvar(self):
        ''' return the next variant of the file, handles headers gracefully ''' 
        line = self.file.readline().strip() 
        while line.startswith('#'):
            # A single # is the column headers
            if line.startswith('##'):
                self.header_buffer.append(line)
            else:
                self.column_header = line
            line = self.file.readline().strip() 
        if line == '':
            return []
        else:
            return Variant(line.split())
 
    # Buffer Full functions
    def max_size(self,stream):
        ''' The buffer contains a maximum number of SNPs '''
        if stream == 'downstream':
            if len(self.buffer) - self.current <= self.max_buffer_size:
                return False
            else:
                return True
        elif stream == 'upstream':
            if self.current > self.max_buffer_size:
                return False
            else:
                return True
            
    def max_window(self,stream):
        ''' The buffer contains a maximum number of SNPs within a window ''' 
        if stream == 'downstream':
            if self.staged != [] and abs(self.staged.pos - self.buffer[self.current].pos) <= self.window_size:
                return False
            else:
                return True
        if stream == 'upstream':
            if abs(self.buffer[0].pos - self.buffer[self.current].pos) > self.window_size:
                return False
            else:
                return True

    def rebuffer(self):
        ''' adjusts buffer so that it satisfies the buffer function ''' 
        if len(self.buffer) == 0: # Beginning condition
            return
        while not self.is_buffer_full('downstream'):
            self.buffer.append(self.staged)
            self.staged = self.nextvar()
        while not self.is_buffer_full('upstream'):
            self.buffer.pop(0)
            self.current -= 1

    def readline(self):
        ''' returns the next variant in the buffer, handles rebuffering '''
        # Beginning condition
        if len(self.buffer) == 0:
            self.buffer.append(self.staged)
            self.staged = self.nextvar()
        elif self.current == len(self.buffer)-1:
            # End condition
            if self.staged == []:
                return self.staged
            self.buffer.append(self.staged)
            self.staged = self.nextvar()
            self.current += 1 # go to next in buffer
        else:        
            self.current += 1 # go to next in buffer
        self.rebuffer()
        return self.buffer[self.current] 

    def header(self):
        return "\n".join(self.header_buffer + [self.column_header])

    def add_header(self,key,val):
        self.header_buffer.append("##{}='{}'".format(key,val))

    def __iter__(self):
        return self

    def __next__(self):
        if self.readline() == []:
            raise StopIteration
        return self.buffer[self.current]


    def __str__(self):
        st = "Buffer:\n-------\n"
        for x,y in enumerate(self.buffer):
            if x == self.current:
                st += str(y)+"<-- Current\n"
            else:
                st += str(y)+"\n"
        st += "Staged:\n-------\n{}".format(self.staged)
        return str(st)
    def __repr__(self):
        return self.__str__()


    def open(self,filename,*args):
        if 'gz' in filename:
            f = gzip.open(filename,'rb')
        else:
            f = open(filename,*args)
        return f
 
