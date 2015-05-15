class VCF(object):
    def __init__(self,vcffile,force=False):
        self.vcffile = open(vcffile,'r')
        self.header = defaultdict(list)
        # keep track of a bunch of indexes
        self.idmap = {}
        self.posmap = defaultdict(dict)
        self.indexmap = []
        # Keep track of samples
        self.samples = []
        # experimental genotype data frame
        self.genotypes = pd.DataFrame()
        # load/create indices
        self.index(force=force)
    def index(self,force=False):
        if not os.path.exists(self.vcffile.name+'.pdx') or force:
            # create the index file
            log('Index file does not exist for {}, indexing now'.format(self.vcffile.name))
            cur_byte = 0
            self.vcffile.seek(0)
            for line in self.vcffile:
                if line.startswith('##'):
                    key,val = line.lstrip('#').rstrip().split('=',1)
                    self.header[key].append(val)
                elif line.startswith("#"):
                    self.samples = line.strip().split()[9:]
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
            pickle.dump((self.idmap,self.posmap,self.samples,self.indexmap),open(self.vcffile.name+".pdx",'wb'))
        else:
            # read the index file
            self.idmap,self.posmap,self.samples,self.indexmap = pickle.load(open(self.vcffile.name+'.pdx','rb')) 

    def iter_variants(self):
        ''' returns variant generator, for iteration. Should be memory efficient '''
        self.vcffile.seek(0)
        return (Variant(line) for line in self.vcffile if not line.startswith('#'))
    def iter_chroms(self):
        return (chrom for chroms in self.posmap.keys())

    def _load_genos(self):
        if not self.genotypes.empty:
            return
        log("Loading genos for {}",self.vcffile.name)
        alleles = []
        ids = []
        for variant in self.iter_variants():  # We only support bi-allelic SNPS right now!!!!
            if not variant.biallelic:
                continue
            alleles.append([variant.chrom,variant.pos] + variant.genos(transform = Allele.vcf2geno))
            ids.append(variant.id)
        self.genotypes = pd.DataFrame(alleles,index=ids,columns=['chrom','pos']+self.samples)
    
    def sample_call_rate(self):
        ''' returns call rate for samples '''
        pass
            
    def sample2index(self,sample):
        ''' return the genotype index for sample '''
        return self.samples.index(sample)

    def __getitem__(self,byte):
        self.vcffile.seek(byte)
        return Variant(self.vcffile.readline())

    def __len__(self):
        return len(self.indexmap)

    @property
    def shape(self):
        return (len(self),len(self.samples))

    @property
    def intervals(self):
        return ["{}:{}".format(chrom,pos) for chrom in self.posmap.keys() for pos in self.posmap[chrom].keys() ]

    def ix(self,index):
        return self[self.indexmap[index]]

    def pos(self,chrom,pos):
        return self[self.posmap[chrom][int(pos)]]

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
            if not variant.biallelic:
                continue
            off,fat,mot = variant.genos([off_i,fat_i,mot_i])
            consistencies.append(Trio(off,fat,mot).consistent())
        if return_raw:
            return consistencies
        return sum(consistencies)/len(consistencies)

    def to_fam(self,filename,fam_id="VCF"):
        with open(filename,'w') as OUT:
            for sample in self.samples:
                print("{}\t{}\t0\t0\t0\t0".format(fam_id,sample),file=OUT) 

    def to_fastTagger(self,prefix=None,maf=0.01,r2=0.99):
        '''outputs to fastTAGGER format. will split by chromosome'''
        if prefix is None:
            prefix = os.path.basename(self.vcffile.name.replace('.vcf',''))
        # Check to see if data is phased
        with open(prefix+'.par','w') as OUT:
            print("\n".join([
               'data_maf={}.maf'.format(prefix),
               'data_matrix={}.matrix'.format(prefix),
               'min_maf={}'.format(maf),
               'window_len=100000',
               'min_r2_1={}'.format(r2),
               'min_r2_2={}'.format(r2),
               'min_r2_3={}'.format(r2),
               'max_len=3',
               'max_merge_window_len=100000',
               'min_bin_size=0',
               'max_covered_times=0',
               'mem_size=0',
               'max_tagSNP_num=0.00',
               'model=MMTagger',
               'output={}_TagSet_MAF_{}_R2_{}'.format(prefix,maf,r2),
            ]),file=OUT)
        with open(prefix+'.maf','w') as MAF, open(prefix+'.matrix','w') as MAT:
            # We need to make two files, data_maf and data_matrix.
            for i,variant in enumerate(self.iter_variants()):
                if i % 50000 == 0:
                    log("processing variant {}",i)
                if not variant.phased:
                    raise Exception("VCF file must be phased: variant {}".format(variant.id))
                if not variant.biallelic:
                    continue
                # figure out minor and major allele freq 
                genos = list(chain.from_iterable([allele.split('|') for allele in variant.genos()]))
                ref_count = genos.count('0')
                alt_count = genos.count('1')
                maf = alt_count/len(genos)
                if maf > 0.5:
                    # MAF is switched around
                    genos = [ '0' if x == '1' else '1' for x in genos ]
                    maf = 1-maf
                    major = variant.alt
                    minor = variant.ref
                else:
                    major = variant.ref
                    minor = variant.alt
                print("{}\t{}\t{}\t{}\t{}".format(variant.id,variant.pos,major,minor,maf),file=MAF)
                print(" ".join(genos),file=MAT)

    def to_hapblock(self,prefix=None):
        if prefix is None:
            prefix = os.path.basename(self.vcffile.name.replace('.vcf',''))
        with open(prefix+'.par','w') as PAR:
            print("\n".join([
                "{}".format(len(self.samples)),
                "{}".format(len(self.indexmap)),
                "100000",
                "1",
                "3 {} {}".format(prefix+'.hapdat',prefix+'.blocks'),
                "2 0.80 0.0499",
                "6 0.80 0.50",
                "2",
                "1 1 1 {}".format(prefix+'.hapids'),
                "1 {}".format(prefix+".happattern"),
                "2"
            ]),file=PAR)
        with open(prefix+'.hapdat','w') as DAT, open(prefix+'.happos','w') as POS, open(prefix+'.hapids','w') as IDS:
            # load the genotype data
            self._load_genos()
            # Sort by position
            self.genotypes = self.genotypes.sort('pos',ascending=True)
            # print location to pos file
            print(len(self.genotypes),"\n".join(map(str,self.genotypes['pos'])),file=POS)
            # print the SNP names 
            print(len(self.genotypes),"\n".join(map(str,self.genotypes.index)),file=IDS)
            print("{}\t{}".format(len(self.samples),len(self.genotypes)))
            genomap = {
                -1 : (0,0), 
                0 : (1,1),
                1 : (1,2),
                2 : (2,2)
            }
            for indiv,geno_series in self.genotypes.transpose()[2:].iterrows():
                print(indiv+"\t",*list(chain(*map(lambda x: genomap[x], geno_series))),sep=" ",file=DAT)

    def to_hmp(self,prefix=None):
        ''' outputs to HapMap format '''
        if prefix is None:
            prefix = self.vcffile.name.replace('.vcf','')
        with open(prefix+'.hmp','w') as HMP:
            # Print Header
            print('\t'.join("rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode".split(' ') + self.samples),file=HMP)
            for i,var in enumerate(self.iter_variants()):
                if i % 50000 == 0:
                    log('Processing variant {}',i)
                if not var.phased:
                    raise Exception("VCF file must be phased: variant {}".format(var.id))
                if not var.biallelic:
                    continue
                # Kaput
                print("\t".join([var.id,"{}/{}".format(var.ref,var.alt),var.chrom,var.pos,'.','-','-','-','-','QC+'] + \
                        [geno.replace('0',var.ref).replace('1',var.alt).replace('|','') for geno in var.genos()]
                ),file=HMP)

    def to_LRTag(self,ld_file,freq_file,prefix=None):
        log("Processing {}",self.vcffile.name)
        if prefix is None:
            prefix = self.vcffile.name.replace('.vcf','')
        # load the pos to id mappings
        idmap = {}
        log("Building idmap")
        for var in self.iter_variants():
            idmap[(str(var.chrom),str(var.pos))] = var.id
        # process ld file
        log('Ouputting lods')
        with open(prefix+'.lod','w') as OUT, open(ld_file,'r') as IN:
            header  = IN.readline()
            for line in IN:
                chrom,pos1,pos2,n,r2 = line.strip().split()
                print("{}\t{}\t{}".format(idmap[chrom,pos1],idmap[chrom,pos2],r2),file=OUT)
        # process freq file
        log('Ouputting mafs')
        with open(prefix+'.maf','w') as OUT, open(freq_file,'r') as IN:
            header = IN.readline()
            for line in IN:
                chrom,pos,n,nchrom,*freqs = line.strip().split()
                if int(n) > 2: # only biallelic
                    continue
                print("{}\t{}".format(self.pos(chrom,pos).id,min(freqs)),file=OUT)


