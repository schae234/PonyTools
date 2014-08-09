#!/usr/bin/env python3

import pandas as pd
import os, errno
import sys
from subprocess import Popen, PIPE
from optparse import OptionParser
from Tools import skipno,log

def main(args):

    parser = OptionParser()
    parser.add_option('--sample-table', dest="sample_table", help="CSV Sample Table")
    parser.add_option('--sample-type', dest="sample_type",   help="CSV Sample Type")
    parser.add_option('--cel-path', dest='cel_path', help="Path of CEL files")
    options,files = parser.parse_args(args) 

    # Read in samples
    samples = pd.read_csv(options.sample_table)
    # Read in sample type
    sample_type = pd.read_table(options.sample_type)
    samples = samples.merge(sample_type[['SampleName','Isolatedfrom']].drop_duplicates(),how="left",left_on="SampleName",right_on="SampleName")
    # Add .CEL to the best array column
    samples['BestArray']+='.CEL'
    
    
    # no samples should have a NULL method oisolation
    assert(any(samples.Isolatedfrom.isnull()) == False)
    
    processes = []
    log('Processing Quality Control')
    for groups,df in samples.groupby(["Design","Isolatedfrom"]):
        # unpack out groups
        design,tissuetype = groups
        log("Found {} samples for group {}",len(df),"/".join(groups))
        # Create the path for this batch
        path = os.path.join(design,tissuetype)
        # make the batch directory
        os.makedirs(path, exist_ok=True)
        # Make the CELs file
        with open(os.path.join(path,"CELS.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            # print out the cels belonging to that batch
            for CEL in df.BestArray.values:
                print(os.path.join(options.cel_path,str(CEL)),file=CELS)
        # EXECUTE!
        stdout = open(os.path.join(path,'QC','qc-geno-stdout.txt'),'w')
        stderr = open(os.path.join(path,'QC','qc-geno-stderr.txt'),'w')
        # Dont run script if the output is already there
        if not os.path.exists(os.path.join(path,'QC',"qc-geno.txt")):
            p = Popen([
                'apt-geno-qc',
                '--analysis-files-path', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/'.format(design,design),
                "--cdf-file",            "Axiom_MNEc2M_{}.v1.cdf".format(design),
                "--chrX-probes",         "Axiom_MNEc2M_{}.v1.chrXprobes".format(design),
                "--chrY-probes",         "Axiom_MNEc2M_{}.v1.chrYprobes".format(design),
                "--target-sketch",       "Axiom_MNEc2M_{}.v1.AxiomGT1.sketch".format(design),
                "--qcc-file",            "Axiom_MNEc2M_{}.v1.qcc".format(design),
                "--qca-file",            "Axiom_MNEc2M_{}.v1.qca".format(design),
                "--female-thresh",       "0.5",
                "--male-thresh",         "0.9",
                '--cel-files',           os.path.join(path,'CELS.txt'),
                '--out-dir',             os.path.join(path,'QC'), 
                '--out-file',            os.path.join(path,'QC','qc-geno.txt'),
                '--log-file',            os.path.join(path,'QC','qc-geno-log.txt')
            ],stdout=stdout,stderr=stderr)
            p.group = "/".join(groups)
            processes.append(p)
    # Wait for apt-geno-qc to finish
    while processes:
        for p in processes[:]:
            if p.poll() is not None:
                print("Process {} ({}) ended with exit code {}".format(p.group,p.pid,p.returncode))
                processes.remove(p)
    ALLDQC = pd.DataFrame()
    # Remove Samples with DQC less than 0.6
    log('Examining DQC and generating priors')
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        design,tissuetype = groups
        path = os.path.join(design,tissuetype)
        # Read in the qc-geno results and filter out DQC below 0.6
        DQCResults = pd.read_table(os.path.join(path,'QC',"qc-geno.txt"),skiprows=skipno(os.path.join(path,'QC',"qc-geno.txt")))
        DQCResults['Design'] = design
        DQCResults['TissueType'] = tissuetype
        ALLDQC = pd.concat([ALLDQC,DQCResults])
        with open(os.path.join(path,"CELS2.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            for i,CEL in DQCResults.iterrows():
                if CEL.axiom_dishqc_DQC < 0.6:
                    print("Dropping {} for {} because {} is < 0.6".format(CEL.cel_files,"/".join(groups),CEL.axiom_dishqc_DQC))
                else:
                    print(os.path.join(options.cel_path,CEL.cel_files),file=CELS)
        # Produce SNP Specific Models for 20K
        path = os.path.join(design,tissuetype,'priors')
        if not os.path.exists(os.path.join(path,"AxiomGT1.snp-posteriors.txt")):
            os.makedirs(os.path.join(path,),exist_ok = True)
            stdout = open(os.path.join(path,'stdout.txt'),'w')
            stderr = open(os.path.join(path,'stderr.txt'),'w')
            p =Popen([
                'apt-probeset-genotype',
                '--analysis', 'brlmm-p-plus.force',
                '--probeset-ids', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallOne/HighPerformingSNPS/{}.txt'.format(design),
                '--write-models',
                "--no-gender-force",
                '--analysis-files-path', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/'.format(design,design),
                "--cdf-file"           ,"Axiom_MNEc2M_{}.v1.cdf".format(design),
                '--set-analysis-name','AxiomGT1',
                '--chip-type', 'Axiom_MNEc2M_{}'.format(design),
                '--summaries',
                "--chrX-probes"        ,"Axiom_MNEc2M_{}.v1.chrXprobes".format(design),
                "--chrY-probes"        ,"Axiom_MNEc2M_{}.v1.chrYprobes".format(design),
                '--cel-files', os.path.join(design,tissuetype,'CELS2.txt'),
                '--out-dir', os.path.join(path), 
                '--log-file', os.path.join(path,'priors-log.txt')
            ],stdout=stdout,stderr=stderr)
            p.group = "/".join(groups)
            processes.append(p)
    # Run model building processes!
    while processes:
        for p in processes[:]:
            if p.poll() is not None:
                print("Process {} (PID {}) ended with exit code {}".format(p.group,p.pid,p.returncode))
                processes.remove(p)
    # Genotype 20K
    log('Genotypeing High Quality SNP set')
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        design,tissuetype = groups
        path = os.path.join(design,tissuetype,'calls-axiom-snp-specific')
        if not os.path.exists(os.path.join(path,"AxiomGT1.calls.txt")):
            os.makedirs(os.path.join(path,),exist_ok = True)
            stdout = open(os.path.join(path,'stdout.txt'),'w')
            stderr = open(os.path.join(path,'stderr.txt'),'w')
            p = Popen([
                'apt-probeset-genotype',
                #'--analysis', 'brlmm-p-plus',
                '--analysis', 'artifact-reduction.ResType=2.Clip=0.4.Close=2.Open=2.Fringe=4.CC=2,quant-norm.target=1000.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.transform=MVA.copyqc=0.00000.wobble=0.05.MS=0.15.copytype=-1.clustertype=2.ocean=0.00001.CSepPen=0.1.CSepThr=4.hints=1.CP=16',
                '--probeset-ids', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallOne/HighPerformingSNPS/{}.txt'.format(design),
                '--analysis-files-path', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/'.format(design,design),
                "--cdf-file"           ,"Axiom_MNEc2M_{}.v1.cdf".format(design),
                '--set-analysis-name','AxiomGT1',
                '--chip-type', 'Axiom_MNEc2M_{}'.format(design),
                '--read-models-brlmmp', os.path.join(design,tissuetype,'priors','AxiomGT1.snp-posteriors.txt'),
                #'--read-models-brlmmp', 'Axiom_MNEc2M_{}.v1.generic_prior.txt'.format(design),
                "--chrX-probes"        ,"Axiom_MNEc2M_{}.v1.chrXprobes".format(design),
                "--chrY-probes"        ,"Axiom_MNEc2M_{}.v1.chrYprobes".format(design),
                '--target-sketch'      ,'Axiom_MNEc2M_{}.v1.AxiomGT1.sketch'.format(design),
                "--set-gender-method"  ,"cn-probe-chrXY-ratio",
                '--qmethod-spec'       ,'med-polish.expon=true',
                "--em-gender"          ,"false",
                "--female-thresh"      ,"0.5",
                "--male-thresh"        ,"0.9", 
                '--cel-files', os.path.join(design,tissuetype,'CELS2.txt'),
                '--out-dir', os.path.join(path), 
                '--log-file', os.path.join(path,'genotype-log.txt')
            ],stdout=stdout,stderr=stderr)
            p.group = "/".join(groups)
            processes.append(p)
    while processes:
        for p in processes[:]:
            if p.poll() is not None:
                print("Process {} (PID {}) ended with exit code {}".format(p.group,p.pid,p.returncode))
                processes.remove(p)
    # Filter Out Samples with low Call Rate
    for groups,df in samples.groupby(['Design','IsolatedFrom']):
        design,tissuetype = groups

    # Filter Out plates with low call Rate

    # Re create Models

    # Re Genotype all SNPs

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
