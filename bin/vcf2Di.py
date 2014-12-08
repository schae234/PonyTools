#!/usr/bin/env python3

import sys
from optparse import OptionParser


def main(args):
	parser = OptionParser()
	parser.add_option('--vcf',help='input vcf file')
	parser.add_option('--out',help='output file prefix')
	parser.add_option('--dummy_qual',default=None, help='Make genotype quality scores up if they are not in the VCF file. specify quality score after flag')
	
	options,args = parser.parse_args(args)

	with open(options.vcf,'r') as IN, open(options.out+'.snpinfo','w') as SNP, open(options.out+'.geno','w') as GEN:
		for line in IN:
			if line.startswith('##'):
				continue
			elif line.startswith('#CHR'):
				# We need the individual ids from the header row
				chrom,pos,id,ref,alt,qual,filter,info,format,*indids = line.strip().split()
				print("Found {} individuals".format(len(indids)),file=sys.stderr)
			else:
				# we are in the genotypes
				chrom,pos,snpid,ref,alt,qual,filter,info,format,*genos = line.strip().split()
				chrom = chrom.replace('chr','')
				# print the snp info
				print("\t".join([snpid,'0',chrom,pos]),file=SNP) 
				# print the genotype information for each individual
				for i,geno in enumerate(genos):
						# gotta extract the genotypes and quality score
						fields = format.split(':')
						geno = geno.split(':')
						allele1,allele2 = geno[fields.index('GT')].replace('|','/').split('/')
						try:
							geno_qual = geno[fields.index('GQ')]
						except ValueError as e:
							if options.dummy_qual:
									geno_qual = str(options.dummy_qual)
							else:
								raise Exception('Genotype Qaulity scores not present in file, see --dummy_qual option.')
						print(','.join([
							snpid,			# snp id
							indids[i],   	# sample id
							'null',			# blank column
						    allele1,				# allele 1
							allele2, 			# allele 2
							'null',
							geno_qual		# genotype quality	
						]),file=GEN)

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
