#!/usr/bin/env python3

from subprocess import call
import sys

vcftools = '/home/grad01/schaefer/bin/vcftools'

def main(argv):
    vcffile,id_file = argv
    # extract the individuals from vcf
    if not os.path.exists(idfile.replace('.txt','')+'.recode.vcf'):
        call([vcftools, 
                '--vcf', vcffile, 
                '--keep', id_file,
                '--recode',
                '--out', id_file.replace('.txt','')]
        )



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
