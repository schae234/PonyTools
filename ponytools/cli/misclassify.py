#!/usr/bin/env python3

import sys
from optparse import OptionParser
from PonyTools.Tools import VCF,log,TagSet
from subprocess import Popen,PIPE
import os

def main(args):
    ''' 
        This script takes in a VCF file (a sample x snp matrux), splits it 
        into --ref and --imp sets. From the --imp set, it further reduces 
        the VCF to a set of --tags. Using the --imp X --tags subset, it uses
        beagle 4 to impute and discovers the miscalssification rate of the 
        
    '''
    parser = OptionParser()
    parser.add_option('--vcf')
    parser.add_option('--ids',default=None)
    parser.add_option('--tags',default=None)
    parser.add_option('--call_thresh',default=0.45,type='float')
    # parse options!
    options,args = parser.parse_args(args)
    processes = []
    # filter imp set into tag set
    log('Generating keep and ref sets' + '#'*50)
    if not os.path.exists('keep_'+os.path.basename(options.ids)+'.recode.vcf'):
        p = Popen([
            'vcftools', 
            '--vcf', options.vcf,
            '--keep', options.ids,
            '--recode',
            '--out', 'keep_' + os.path.basename(options.ids)
        ])
        p.desc = 'vcftools: keep'
        processes.append(p)
    if not os.path.exists('ref_'+os.path.basename(options.ids)+'.recode.vcf'):
        p = Popen([
            'vcftools', 
            '--vcf', options.vcf,
            '--remove', options.ids,
            '--recode',
            '--out', 'ref_' + os.path.basename(options.ids)
        ])
        p.desc = 'vcftools: remove' 
        processes.append(p)
    # Wait for vcftools to be done
    wait(processes)
    # filter to only tags for keep set
    if not os.path.exists('tags_'+os.path.basename(options.ids)+'.recode.vcf'):
        log('Genertaing Tag VCF'+ '#'*50)
        p = Popen([
            'vcftools',
            '--vcf', 'keep_' + os.path.basename(options.ids) + '.recode.vcf',
            '--snps', options.tags,
            '--recode',
            '--out', 'tags_' + os.path.basename(options.ids)
        ])
        p.desc = 'vcftools: tags'
        processes.append(p)
        wait(processes)
    # phase the reference pops
    log('Phasing reference set' + '#'*50)
    if not os.path.exists('ref_phased_'+os.path.basename(options.ids)+'.recode.vcf'):
        p = Popen([
            'java', 
            '-jar', '/project/mccuelab/shared/bin/b4.r1274.jar',
            'gt=ref_'+os.path.basename(options.ids)+'.recode.vcf',
            'out=ref_phased_'+os.path.basename(options.ids)+'.recode'
        ])
        p.desc = 'ref phasing'
        processes.append(p)
        wait(processes)
        unzip('ref_phased_'+os.path.basename(options.ids)+'.recode.vcf.gz',processes)
    log('Imputing Tags' + '#'*50)
    if not os.path.exists('imputed_'+os.path.basename(options.ids)+'.recode.vcf'):
        p = Popen([
            'java', 
            '-jar', '/project/mccuelab/shared/bin/b4.r1274.jar',
            'gt=tags_'+os.path.basename(options.ids)+'.recode.vcf',
            'ref=ref_phased_'+os.path.basename(options.ids)+'.recode.vcf',
            'out=imputed_'+os.path.basename(options.ids)+'.recode',
        ])
        p.desc = 'tag imputation'
        processes.append(p)
        wait(processes)
        unzip('imputed_'+os.path.basename(options.ids)+'.recode.vcf.gz',processes)

    # Now do what you came here to do!
    tags = TagSet('tags_'+os.path.basename(options.ids)+'.recode.vcf')
    print(os.path.basename(options.ids),
        tags.misclassification(
            VCF("imputed_"+os.path.basename(options.ids)+".recode.vcf"),
            VCF("keep_"+os.path.basename(options.ids)+".recode.vcf"),
            call_thresh=options.call_thresh
        )
    )


def wait(processes):
    '''waits on a list of processes'''  
    while processes:
        for p in processes:
            if p.poll() is not None:
                print("Process {} ({}) ended with exit code {}".format(p.desc,p.pid,p.returncode))
                processes.remove(p) 

def unzip(file,processes):
    p = Popen(['gunzip','-f',file])    
    p.desc = 'unzipping'
    processes.append(p)
    wait(processes)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
