#!/usr/bin/env python3

import os
import sys
import numpy
import re
import bisect
import time
from optparse import OptionParser,OptionGroup
from collections import defaultdict
from MNEcTools.VCFBuffer import VCFBuffer


class VCFFilter(object):
    def __init__(self):
        self.percentile = None
        self.field = None
        self.qual = None
        self.thresh = None
        self.union_positions = None
        self.map_positions = None
        self.exclude_map = None
        self.out = sys.stdout
        self.log_file = sys.stderr
        self.skip_chroms = None
        self.keep_chroms = None
        self.num_call_sets = None
        self.repeats = None
        self.allelic = None
        self.emit_type = "all"
        self.fasta = None
        self.priority = '0'
        # VCF Position paramaters
        self.buffer_function = "max_size"
        self.buffer_size = 10
        self.within = None
        self.not_within = None
        self.max_ld = None
        self.ld = None
        self.preference_set = None
        self.preference_field = None
        # Set counters
        self.num_within = 0
        self.num_allelic = 0
        self.num_repeat = 0
        self.num_field = 0
        self.num_ld = 0
        self.num_pref_set = 0
        self.num_map = 0
        self.buffer = None
        self.start_time = time.ctime() 


    
    def process(self,vcf_file):
        ''' This method performs the filtering. It assumes you set up the filters ahead of time '''
        self.log("Processing: {}".format(vcf_file))
        # Filter the vcf file
        self.Varbuffer = VCFBuffer(vcf_file,self.buffer_function,self.buffer_size)
        # Print the varaint buffer   
        self.Varbuffer.add_header('vcf_qual_percentile',' '.join(sys.argv))
        self.emit_header()
        for ln,variant in enumerate(self.Varbuffer):
            if ln % 100000 == 0:
                self.log("\tprocessing variant: {}".format(ln))
                self.filter_totals()
                self.out.flush()
            self.filter_by_property(variant)
            self.filter_by_proximity(self.Varbuffer)
            self.emit(variant)
                
    def filter_by_proximity(self,Varbuffer):           
        variant = Varbuffer.buffer[Varbuffer.current]
        if self.not_within:
            if Varbuffer.current == 0:
                upstream = float('Inf')
            else:
                upstream = Varbuffer.buffer[Varbuffer.current-1].pos
            if (len(Varbuffer.buffer) == Varbuffer.current+1
                or Varbuffer.buffer[Varbuffer.current+1] == []): 
                downstream = float('Inf')
            else:
                downstream = Varbuffer.buffer[Varbuffer.current+1].pos
            if ( abs(variant.pos - upstream )<= self.not_within 
              or abs(variant.pos - downstream)<= self.not_within):
                variant.add_filter("WITHIN")
                self.num_within += 1
        if self.ld:
            if variant.pos in self.ld[variant.chrom] and "MAX_LD" not in variant.filters:
                # if in LD, we want to filter further
                # Get info about variants in LD
                in_ld = self.ld[variant.chrom][variant.pos].keys()
                buff_pos = map(lambda x: x.pos, self.Varbuffer.buffer)
                ld_ind = [ buff_pos.index(x) for x in buff_pos if x in in_ld ]
                ld_field = [ self.Varbuffer.buffer[x].info[self.preference_field] for x in ld_ind ]
                # collapse down all the sets of the variants in ld
                ld_sets = set([item for item in self.Varbuffer.buffer[x].info['set'].split('-') for x in ld_ind ])
                # Get info about the variant in question
                variant_field = variant.info[self.preference_field]
                variant_sets = set(variant.info['set'].split('-'))
                # If the AF == 1.0 and also in the thoroughbreds, its probably a FP
                if float(variant.info['AF']) == 1.0 and ('variant2' in variant_sets or 'variant6' in variant_sets):
                    variant.add_filter('REFERENCE_MISCALL')
                    return
                    # Make preference of certain sets (ponies or drafts)
                if self.preference_set:
                    if len(variant_sets.intersection(self.preference_set)) > len(ld_sets.intersection(self.preference_set)):
                        # mark all ld variants as filtered
                        [ self.Varbuffer.buffer[x].add_filter("MAX_LD_PREF_SET") for x in ld_ind ]
                        self.num_pref_set += len(ld_ind)
                        return
                    elif len(variant_sets.intersection(self.preference_set)) < len(ld_sets.intersection(self.preference_set)):
                        # mark the variant in question as filtered
                        variant.add_filter("MAX_LD_PREF_SET")
                        self.num_pref_set += 1
                        return
                # pref sets are seen equally in both
                # prioritize by ld
                if all([variant_field > x for x in ld_field]):
                    # Add filter to all ld variants
                    [ self.Varbuffer.buffer[x].add_filter("MAX_LD") for x in ld_ind ]
                    self.num_ld += len(ld_ind)
                else:
                    variant.add_filter("MAX_LD")
                    self.num_ld += 1

     
    def filter_by_property(self,variant):
        ''' Filter variant due to certain properties '''
        if self.allelic and len(variant.alt.split(',')) >= self.allelic:
            variant.add_filter("ALLELIC")
            self.num_allelic += 1
            return
        if self.skip_chroms and variant.chrom in self.skip_chroms:
            variant.add_filter("SKIP_CHROM")
        if self.keep_chroms and variant.chrom not in self.keep_chroms:
            variant.add_filter("NOT_KEEP_CHROM")
        if self.qual and variant.qual < self.qual:
            variant.add_filter("LOW_QUAL_THRESH")
        if self.field and self.thresh:
            for i,field in enumerate(self.field):
                try:
                    if float(variant.info[field]) < self.thresh[i]:
                        variant.add_filter("LOW_"+field.upper())
                        self.num_field += 1
                except Exception as e: 
                    variant.add_filter("NO_FIELD"+field.upper())
                    #self.log("No {} field for variant {} {}".format(
                    #    field,variant.chrom,variant.pos
                    #))
        if self.repeats and self.in_repeat_region(variant.chrom,variant.pos):
            variant.add_filter("REPEAT")
            self.num_repeat += 1
        if self.union_positions and variant.pos not in self.union_positions[variant.chrom]:
            variant.add_filter("NON_INTERSECT")
        if self.map_positions:
            if self.exclude_map and variant.pos in self.map_positions[variant.chrom]:
                variant.add_filter("MAP")
                self.num_map += 1
            elif not self.exclude_map and variant.pos not in self.map_positions[variant.chrom]:
                variant.add_filter("MAP")
                self.num_map += 1
        if self.num_call_sets:
            try:
                sets = variant.info['set'].split('-')
                if len(sets) < self.num_call_sets or 'FilteredInAll' in sets:
                    variant.add_filter("LOW_CALL_SET")
            except Exception as e:
                self.log("No set found for {} {}".format(fields[0],fields[1]))
        if variant.filters == '.':
            variant.filters = "PASS"


    def populate_intersect(self,vcf_file):
        self.log("Reading in union file: {}".format(vcf_file))
        positions = defaultdict(list)
        with open(vcf_file,'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                line = line.rstrip()
                line = line.split()
                positions[line[0]].append(int(line[1]))
        num_vars_seen = 0
        for chr in positions.keys():
            # covert to set for fast intersection checks
            positions[chr] = set(positions[chr])
            num_vars_seen += len(positions[chr])
        self.log("Found {} variants in union file".format(num_vars_seen))
        self.union_positions = positions

    def get_percentile_cutoff(self,vcf_file,perc):
        ''' Return the field cutoff for a specific field cutoff ''' 
        with open(vcf_file,'r') as file:
            quals = numpy.array(
                [float(line.split()[5]) for line in file if not line.startswith("#") ]
            )
        qual_cutoff = numpy.percentile(quals,perc)
        return qual_cutoff 


    def populate_map(self,map_file):
        positions = defaultdict(list)
        with open(map_file,'r') as file:
            for line in file:
                line=line.strip()
                chr,id,cm,pos = line.split()
                positions['chr'+chr].append(int(pos))
        for chr in positions.keys():
            # covert to set for fast intersection checks
            positions[chr] = set(positions[chr])
        self.map_positions = positions
    
    def populate_repeat(self,repeat_file):
        self.log("Reading in repeat file: {}".format(repeat_file))
        ranges = defaultdict(lambda: defaultdict(list))
        with open(repeat_file,'r') as file:
            for line in file:
                chr,start,end = line.strip().split() 
                ranges[chr]['start'].append(int(start))
                ranges[chr]['end'].append(int(end))
        self.repeats = ranges
    
    def in_repeat_region(self,chrom,position):
        if chrom not in self.repeats:
            return False # Be on the safe side
        index = min(bisect.bisect_left(self.repeats[chrom]['start'],position)-1,0)
        if self.repeats[chrom]['start'][index] < position and self.repeats[chrom]['end'][index] > position:
            return True
        else:
            return False

    def populate_ld(self,ld_file):
        '''Processes pairwise ld file produced by vcftools'''
        self.log("Reading in LD File")
        self.ld = defaultdict(lambda:defaultdict(lambda: defaultdict(int)))
        cur_chrom = None
        with open(ld_file,'r') as file:
            header = file.readline()
            for line in file:
                chr,pos1,pos2,n_ind,r2 = line.strip().split()
                if float(r2) >= self.max_ld and int(n_ind) > 3:
                    # we only store the top half of the matrix
                    self.ld[chr][int(pos1)][int(pos2)] = (int(n_ind),float(r2))
                    #self.ld[chr][int(pos2)][int(pos1)] = (int(n_ind),float(r2))
                if chr != cur_chrom:
                    self.log("\tprocessing chromosome {}".format(chr))
                    cur_chrom = chr

    # IO Methods
    def log(self,*args):
        print(time.ctime(),' - ',*args,file=self.log_file)
        self.log_file.flush()

    def emit_header(self):
        # skip header for some emit types
        if self.emit_type not in  ['all','pass','filtered']:
            pass
        else:
            print(self.Varbuffer.header(),file=self.out) 

    def emit(self,variant):
        # emit fields
        if self.emit_type == 'info':
            print(*self.info(variant),file=self.out)
        # emit stats about each variant
        elif self.emit_type == 'stats':
            print(*self.stats(variant),file=self.out) 
        # emit only passing variants
        elif self.emit_type == 'map':
            # Map format follows: chrom id cm pos
            print("{}\t{}\t{}\t{}".format(variant.chrom,variant.id,'0',variant.pos), file=self.out)
        elif self.emit_type == 'snpchip':
            print("\t".join([
                    "EquCab",                               # Organism
                    variant.MNEcId,                         # SNPID 
                    self.flank_seq(variant,flank_len=35),   # SEQ
                    self.priority,                                    # SNP_Priority
                    variant.chrom,                          # CHR
                    str(variant.pos),                            # POS
                    variant.chr_type,                       # CHR_TYPE
                    '0'                                     # SNP_VAL (validations) 
                ]), 
            file=self.out)
        elif self.emit_type == 'pass' and (len(variant.filters) > 1 or 'PASS' not in variant.filters):
            return None
        # emit only filtered variants
        elif self.emit_type == 'filtered' and (fields == 'PASS' or fields[6] == '.'):
            return None
        # base case
        else:
            print(variant,file=self.out)

    def stats(self,variant):
        return([
            variant.info['VQSLOD'],
            variant.qual,
            variant.info['AF'],
            abs(self.Varbuffer.buffer[self.Varbuffer.nexti].pos - self.Varbuffer.buffer[self.Varbuffer.previ].pos),
            'variant' in variant.info['set'].split('-'),
            'variant9' in variant.info['set'].split('-'),
        ])
    def info(self,variant):
        return([
            variant.MNEcId,
            variant.chrom,
            str(variant.pos),
            "T" if any([pref_set in variant.info['set'] for pref_set in self.preference_set ]) else "F"
        ])

    def filter_totals(self):
        self.log(''' 
            Filter totals:
            within {}
            allelic {}
            repeat {}
            field {}
            ld {}
            preference set {}
            num map {}
        '''.format(self.num_within,self.num_allelic,self.num_repeat,self.num_field,self.num_ld,self.num_pref_set,self.num_map)) 


    #------------ SNP CHIP STUFF -----------------------#
    def process_fasta(self, fasta_file):
        self.log("Reading in Fasta")
        genome = defaultdict(list)
        current_chrom = 0
        with open(fasta_file,'r') as fasta:
            for line in fasta:
                line = line.rstrip()
                if line.startswith(">"):
                    current_chrom = line.split()[0].replace(">","")
                    if self.verbose:
                        self.log("\tOn chromosome {}".format(current_chrom))
                else:
                    # speed things up
                    genome[current_chrom].append(line.upper())
        self.fasta = {
            key:"".join(val) for key,val in genome.items()
        }
    
    def flank_seq(self,variant,flank_len=80):
        assert self.fasta 
        index = variant.pos-1
        if self.fasta[variant.chrom][index] != variant.ref: # and fasta[chr][ind] != alt:
            exit("SNP didnt match reference for variant {}".format(variant.MNEcId))
        flank = self.fasta[variant.chrom][index-flank_len:index] \
            + "[" + "/".join(sorted([variant.ref,variant.alt])) + "]" \
            + self.fasta[variant.chrom][index+1:index+flank_len+1]
        return flank



def main(argv):
    parser = OptionParser(usage="Usage: %prog [options] vcf1 [...vcfn]",
        description="This program filters VCF files according to various number of options."
    )
    # General Options/Filters
    gen_filters = OptionGroup(parser,"General Options/Filters")
    gen_filters.add_option(
        "--filter",default=None, type=str,
            help="filter out variants that do not match based on FILTER field, e.g. LowQual")
    gen_filters.add_option(
        "--field", default=None, type=str, action="append",
            help="Used with percentile or threshold arguments to filter out based on a FIELD arg")
    gen_filters.add_option(
        "--qual", default=None,  type=int,
            help="Used with percentile or threshold arguments to filter out based on a QUAL arg")
    gen_filters.add_option(
        "--thresh", default=None,type=float,action="append",
            help='filter out --field values under this threshold')
    gen_filters.add_option(
        "--percentile", default = None, type = int,
            help="filter out --field values under this percentile")
    gen_filters.add_option(
        '--skip_chrom',action="append",type=str,
            help="Skip variants on chroms. Can be specified multiple times")
    gen_filters.add_option(
        '--keep_chrom',action="append",type=str,
            help="keep variants on chroms. Can be specified multiple times")
    gen_filters.add_option(
        '--num_call_sets',type=float,
            help="the minimum number of call sets a variant must be called in")
    gen_filters.add_option(
        '--allelic',type=int,default=None,
            help="Filter out variants which have more alleles than X")
    parser.add_option_group(gen_filters)
    # Options that are are proximal
    prox_filter = OptionGroup(parser,"Proximity Filters","Filters which depend on the values of other variants.")
    prox_filter.add_option(
        '--buffer_function',type=str,default=None,
            help="Determines how program will buffer proximal variants, choose max_size or max_window")
    prox_filter.add_option(
        '--buffer_size',type=int,default=None,
            help="Comparison window in base pairs.")
    prox_filter.add_option(
        '--max_ld',type=float,default=None,
            help="Max linkage disequilibrium allowed within a window.")
    prox_filter.add_option(
        '--not_within',type=int,default=None,
            help="the variant is not within X base pairs to another variant, if so, both get filtered")
    prox_filter.add_option(
        '--preference_set',type=str,default=None,action='append',
            help="If any tie breaking needs to happen, take preference with this set. May specify multiple sets")
    prox_filter.add_option(
        '--preference_field',type=str,default=None,
            help='If any tie breaking needs to happend, take preference to the variant with the higher field. Specify once.')
    parser.add_option_group(prox_filter)
    # Options that require extraneous files
    extra_filters = OptionGroup(parser,"Extraneous Filters","Filters which require extraneous files.")
    extra_filters.add_option(
        '--repetitive',type=str,
            help="Filters out variants in repetitive regions. Supply a three columned file with chr,start,end regions")
    extra_filters.add_option(
        '--map_file',type=str,default=None,
            help="Only include snps in Plink map file")
    extra_filters.add_option(
        '--exclude_map_file',type=str,default=None,
            help='Exclude SNPs found in map file')
    extra_filters.add_option(
        "--intersect_vcf", default = None, type=str, 
            help="Only match variants within another vcf")
    extra_filters.add_option(
        "--exclude_vcf",default=None, type=str,
            help="Only match variants not within another vcf")
    extra_filters.add_option(
        "--ld_file",default = None, type=str, 
            help="A file containing LD measurements produced by vcftools, 5 cols: chr, pos1, pos2, #indivs, r^2")
    extra_filters.add_option(
        '--fasta',default=None,type=str,
            help="a fasta file containing the sequence for the genome"
        )
    # Options related to output and formatting
    extra_filters.add_option(
        "--verbose", default = False, action="store_true",
            help="how much do you want to know?")
    extra_filters.add_option(
        '--emit',type=str,default="all",
            help="Controls which variants to print out. One of all, pass, or filtered or stats")
    extra_filters.add_option(
        '--out', default=sys.stdout,
             help="Specify the output file")
    extra_filters.add_option(
        '--priority', default='0', type=str,
        help='priority for variants in output files'
    )
    parser.add_option_group(extra_filters)

    options, args = parser.parse_args(argv) 

    # Create an empty filter and add based on options
    flt = VCFFilter()
    flt.emit_type = options.emit
    flt.verbose = options.verbose
    if options.buffer_function:
        if not options.buffer_size:
            exit('--buffer_function requires --buffer_size')
        flt.buffer_function = options.buffer_function
        flt.buffer_size = options.buffer_size
    # Check to see if we are performing proximity checks
    flt.log("Starting Analysis...")
    # Check to see if we are performing union checks
    if options.intersect_vcf != None:
        flt.populate_intersect(options.intersect_vcf)
    # Check for map based filters
    if options.map_file != None:
        flt.populate_map(options.map_file)
    if options.exclude_map_file != None:
        flt.populate_map(options.exclude_map_file)
        flt.exclude_map = True
    # Check for chrom based filters
    if options.skip_chrom and options.keep_chrom:
        exit("Cannot specify --keep_chrom AND --skip_chrom")
    if options.skip_chrom:
        flt.skip_chroms = set(options.skip_chrom)
    if options.keep_chrom:
        flt.keep_chroms = set(options.keep_chrom)
    # Check for field based filters, can be either a thresh OR percentile, not both
    if options.field and (options.field or options.percentile):
        flt.field = options.field
        flt.thresh = options.thresh
        flt.percentile = options.percentile
    # Check the threshold based filter
    if options.qual:
        if not options.thresh and not options.percentile:
            exit("--qual requires --thresh or --percentile")
        flt.qual = options.qual
        flt.thresh = options.thresh
    # Check for call set based filters
    if options.num_call_sets:
        flt.num_call_sets = options.num_call_sets
    # Check for percentile based filters
    if options.percentile:
        flt.qual = "QUAL"
        flt.thresh = flt.get_percentile_cutoff(args[0],options.percentile)
        print("QUAL THRESH: {}".format(flt.thresh),file=sys.stderr)
    if options.repetitive:
        flt.populate_repeat(options.repetitive)
    if options.out != sys.stdout:
        flt.out = open(options.out,'w')
    # Check for FILTER based filters
    if options.filter:
        flt.filter = options.filter
    if options.allelic:
        flt.allelic = options.allelic
    # Check for LD based filters
    if options.ld_file:
        if not options.max_ld:
            exit("--ld_file option requires --max_ld parameter")
        flt.max_ld = options.max_ld
        flt.populate_ld(options.ld_file)
    if options.preference_set:
            flt.preference_set = set(options.preference_set)
    if options.preference_field:
        flt.preference_field = options.preference_field
    if options.fasta:
        flt.process_fasta(options.fasta)
    if options.priority:
        flt.priority = options.priority

    # Process input files
    for filename in args:
        # Check if we are filtering proximity:
        if options.not_within:
            flt.not_within = options.not_within
        flt.process(filename)
    if flt.verbose:
        flt.filter_totals()
    flt.log("... analysis ended")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))





