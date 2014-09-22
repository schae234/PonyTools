#!/usr/bin/env python3

import pandas as pd
import os, errno
import sys
import numpy as np
from subprocess import Popen, PIPE
from optparse import OptionParser
from Tools import skipno,log

def main(args):

    parser = OptionParser()
    parser.add_option('--sample-table', dest="sample_table",    help="CSV Sample Table")
    parser.add_option('--sample-type',  dest="sample_type",     help="CSV Sample Type")
    parser.add_option('--cel-path',     dest='cel_path',        help="Path of CEL files")
    parser.add_option('--analysis-path',dest='ana_path',        help="Path of Analysis Files (.cdf, .{X,Y}probes, *.qcc, etc.)")
    options,files = parser.parse_args(args) 

    # Read in samples
    samples = pd.read_csv(options.sample_table)
    # Read in sample type
    sample_type = pd.read_table(options.sample_type)
    samples = samples.merge(sample_type[['SampleName','Isolatedfrom']].drop_duplicates(),how="left",left_on="SampleName",right_on="SampleName")
    # Add .CEL to the best array column
    samples['BestArray']+='.CEL'
    samples.index = samples['BestArray']
    samples['DQC'] = np.nan
    samples['20K_CR'] = np.nan
    samples['Geno_CR'] = np.nan
    samples['DropReason'] = ''
    
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
        os.makedirs(os.path.join(path,'QC'),exist_ok=True)
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
                '--analysis-files-path', options.ana_path,
                "--cdf-file",            "Axiom_MNEc2M_{}.r1.cdf".format(design),
                "--chrX-probes",         "Axiom_MNEc2M_{}.r1.chrXprobes".format(design),
                "--chrY-probes",         "Axiom_MNEc2M_{}.r1.chrYprobes".format(design),
                "--target-sketch",       "Axiom_MNEc2M_{}.r1.AxiomGT1.sketch".format(design),
                "--qcc-file",            "Axiom_MNEc2M_{}.r1.qcc".format(design),
                "--qca-file",            "Axiom_MNEc2M_{}.r1.qca".format(design),
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
        log("Found {} samples for group {}",len(df),"/".join(groups))
        path = os.path.join(design,tissuetype)
        # Read in the qc-geno results and filter out DQC below 0.6
        DQCResults = pd.read_table(os.path.join(path,'QC',"qc-geno.txt"),skiprows=skipno(os.path.join(path,'QC',"qc-geno.txt")))
        ALLDQC = pd.concat([ALLDQC,DQCResults])
        with open(os.path.join(path,"CELS2.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            for i,CEL in DQCResults.iterrows():
                if CEL.axiom_dishqc_DQC < 0.6:
                    samples.loc[CEL.cel_files,'DropReason'] = 'DQC'            
                    print("Dropping {} for {} because {} is < 0.6".format(CEL.cel_files,"/".join(groups),CEL.axiom_dishqc_DQC))
                else:
                    print(os.path.join(options.cel_path,CEL.cel_files),file=CELS)
                samples.loc[CEL.cel_files,'DQC'] = CEL.axiom_dishqc_DQC           
                
        # Produce SNP Specific Models for 20K
        # Genotype 20K
        path = os.path.join(design,tissuetype,'20K')
        if not os.path.exists(os.path.join(path,"AxiomGT1.report.txt")):
            os.makedirs(os.path.join(path,),exist_ok = True)
            stdout = open(os.path.join(path,'stdout.txt'),'w')
            stderr = open(os.path.join(path,'stderr.txt'),'w')
            if len(df) > 96:
                priors = "Axiom_MNEc2M_{}.r1.generic_prior.txt".format(design)
            else:
                log("{} has less than 96 individuals, using snp specific priors",groups)
                priors = "Axiom_MNEc2M_{}.r1.AxiomGT1.models".format(design)
            p =Popen([
                'apt-probeset-genotype',
                '--analysis-files-path', options.ana_path,
                #From XML File Step 1
                "--set-analysis-name"  ,"AxiomGT1" ,
                "--chip-type"          ,"Axiom_MNEc2M_{}".format(design) ,
                "--chip-type"          ,"Axiom_MNEc2M_{}.r1".format(design) ,
                "--analysis"           ,"artifact-reduction.ResType=2.Clip=0.4.Close=2.Open=2.Fringe=4.CC=2,quant-norm.target=1000.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.transform=MVA.copyqc=0.00000.wobble=0.05.MS=0.15.copytype=-1.clustertype=2.ocean=0.00001.CSepPen=0.1.CSepThr=4.hints=1.CP=16",
                "--qmethod-spec"       ,"med-polish.expon=true" ,
                "--read-models-brlmmp" , priors,
                # These probes are the high quality probes ~20,000
                "--probeset-ids"       ,"Axiom_MNEc2M_{}.r1.step1.ps".format(design) ,
                "--cdf-file"           ,"Axiom_MNEc2M_{}.r1.cdf".format(design) ,
                "--chrX-probes"        ,"Axiom_MNEc2M_{}.r1.chrXprobes".format(design) ,
                "--chrY-probes"        ,"Axiom_MNEc2M_{}.r1.chrYprobes".format(design) ,
                "--target-sketch"      ,"Axiom_MNEc2M_{}.r1.AxiomGT1.sketch".format(design) ,
                "--no-gender-force"    ,"true" ,
                "--set-gender-method"  ,"cn-probe-chrXY-ratio" ,
                "--em-gender"          ,"false" ,
                "--female-thresh"      ,"0.5" ,
                "--male-thresh"        ,"0.9" ,
                # End XML File
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
    # Filter Out Samples with low Call Rate
    log('Filtering on 20K call rate and Producing all genotypes')
    # Re Genotype all SNPs
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        design,tissuetype = groups
        log("Found {} samples for group {}",len(df),"/".join(groups))
        path = os.path.join(design,tissuetype)
        # Read in the qc-geno results and filter out DQC below 0.6
        CRResults = pd.read_table(os.path.join(path,'20K',"AxiomGT1.report.txt"),skiprows=skipno(os.path.join(path,'20K',"AxiomGT1.report.txt")))
        #ALLDQC = pd.concat([ALLDQC,DQCResults])
        with open(os.path.join(path,"CELS3.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            for i,CEL in CRResults.iterrows():
                if CEL.call_rate < 97:
                    samples.loc[CEL.cel_files,'DropReason'] = '20K_CR'            
                    print("Dropping {} for {} because {} is < 0.97".format(CEL.cel_files,"/".join(groups),CEL.call_rate))
                else:
                    print(os.path.join(options.cel_path,CEL.cel_files),file=CELS)
                samples.loc[CEL.cel_files,'20K_CR'] = CEL.call_rate            
        # Genotype ALL
        path = os.path.join(design,tissuetype,'geno')
        if not os.path.exists(os.path.join(path,"AxiomGT1.report.txt")):
            os.makedirs(os.path.join(path,),exist_ok = True)
            stdout = open(os.path.join(path,'stdout.txt'),'w')
            stderr = open(os.path.join(path,'stderr.txt'),'w')
            if len(df) > 96:
                priors = "Axiom_MNEc2M_{}.r1.generic_prior.txt".format(design)
            else:
                log("{} has less than 96 individuals, using snp specific priors",groups)
                priors = "Axiom_MNEc2M_{}.r1.AxiomGT1.models".format(design)
            p = Popen([
                'apt-probeset-genotype',
                '--analysis-files-path', options.ana_path,
                #From XML File Step 1
                "--set-analysis-name"  ,"AxiomGT1" ,
                "--chip-type"          ,"Axiom_MNEc2M_{}".format(design) ,
                "--chip-type"          ,"Axiom_MNEc2M_{}.r1".format(design) ,
                "--analysis"           ,"artifact-reduction.ResType=2.Clip=0.4.Close=2.Open=2.Fringe=4.CC=2,quant-norm.target=1000.sketch=50000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.transform=MVA.copyqc=0.00000.wobble=0.05.MS=0.15.copytype=-1.clustertype=2.ocean=0.00001.CSepPen=0.1.CSepThr=4.hints=1.CP=16",
                "--qmethod-spec"       ,"med-polish.expon=true" ,
                "--read-models-brlmmp" , priors,
                "--cdf-file"           ,"Axiom_MNEc2M_{}.r1.cdf".format(design) ,
                "--chrX-probes"        ,"Axiom_MNEc2M_{}.r1.chrXprobes".format(design) ,
                "--chrY-probes"        ,"Axiom_MNEc2M_{}.r1.chrYprobes".format(design) ,
                "--target-sketch"      ,"Axiom_MNEc2M_{}.r1.AxiomGT1.sketch".format(design) ,
                "--no-gender-force"    ,"true" ,
                "--set-gender-method"  ,"cn-probe-chrXY-ratio" ,
                "--em-gender"          ,"false" ,
                "--female-thresh"      ,"0.5" ,
                "--male-thresh"        ,"0.9" ,
                "--write-models",
                "--summaries",
                # End XML File
                '--cel-files', os.path.join(design,tissuetype,'CELS3.txt'),
                '--out-dir', os.path.join(path), 
                '--log-file', os.path.join(path,'geno-log.txt')
            ],stdout=stdout,stderr=stderr)
            p.group = "/".join(groups)
            processes.append(p)
    # Run model building processes!
    while processes:
        for p in processes[:]:
            if p.poll() is not None:
                print("Process {} (PID {}) ended with exit code {}".format(p.group,p.pid,p.returncode))
                processes.remove(p)
    log('Adding Geno call rate')
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        design,tissuetype = groups
        path = os.path.join(design,tissuetype)
        # Read in the qc-geno results and filter out DQC below 0.6
        CRResults = pd.read_table(os.path.join(path,'geno',"AxiomGT1.report.txt"),skiprows=skipno(os.path.join(path,'geno',"AxiomGT1.report.txt")))
        #ALLDQC = pd.concat([ALLDQC,DQCResults])
        with open(os.path.join(path,"CELS4.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            for i,CEL in CRResults.iterrows():
                if CEL.call_rate < 97:
                    samples.loc[CEL.cel_files,'DropReason'] = 'Geno_CR'            
                    print("Dropping {} for {} because {} is < 0.97".format(CEL.cel_files,"/".join(groups),CEL.call_rate))
                else:
                    print(os.path.join(options.cel_path,CEL.cel_files),file=CELS)
                samples.loc[CEL.cel_files,'Geno_CR'] = CEL.call_rate            
        # Produce SNP Specific Models for 20K
        # Genotype 20K
        path = os.path.join(design,tissuetype,'geno')
    samples.to_csv('Samples.csv',index=None)
    # Filter Out plates with low call Rate
    for groups,df in samples.groupby(['Design','Isolatedfrom','RunPlateName']):
        design,ifrom,plate = groups
        num_blood = sum([x == 'blood' for x in df.Isolatedfrom])
        num_hair = sum([x == 'hair' for x in df.Isolatedfrom])
        plate_pass_rate = (sum([ x not in ['20K_CR','DQC'] for x in df.DropReason])/len(df))*100
        mean_call_rate = df[[x not in ['20K_CR','DQC'] for x in df.DropReason]]['20K_CR'].mean()
        print("Analyzing {}".format(plate))
        print("\t{} blood".format(num_blood))
        print("\t{} hair".format(num_hair))
        if (plate_pass_rate < 95 and ifrom =='blood') or (plate_pass_rate < 93 and ifrom =='hair'):
            print("\tplate: {} has poor pass rate: {}".format(plate,plate_pass_rate))
        if mean_call_rate < 99:
            print("\tplate: {} has poor call rate: {}".format(plate,mean_call_rate))
    # Run R script
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        design,tissuetype = groups
        path = os.path.join(design,tissuetype,'geno')
        if not os.path.exists(os.path.join(path,'Ps.performance.txt')):
            p = Popen([
                os.path.expanduser('~/Documents/Codes/MNEcTools/Metrics.R'),
                os.path.join(path,"AxiomGT1.snp-posteriors.txt"),
                os.path.join(path,"AxiomGT1.calls.txt"),
                os.path.join(path,"AxiomGT1.metrics.txt"),
                os.path.join(path,"AxiomGT1.metrics.txt"),
                os.path.join("/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/R1/Axiom_MNEc2M_Analysis.r1/Axiom_MNEc2M_{}.r1.ps2snp_map.ps".format(design)),
                os.path.join(path),

            ],stdout=stdout,stderr=stderr)
            p.group = '/'.join(groups)
            processes.append(p)
    while processes:
        for p in processes[:]: 
            if p.poll() is not None:
                print("Process {} (PID {}) ended with exit code {}".format(p.group,p.pid,p.returncode))
                processes.remove(p)

    allMetrics = pd.DataFrame()
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        log("processing {}","/".join(groups))
        design,tissuetype = groups
        path = os.path.join(design,tissuetype,'geno')
        metrics = pd.read_table(os.path.join(path,'Ps.performance.txt'))
        metrics['design'] = design
        metrics['tissue'] = tissuetype
        snp_info = pd.read_table(
            "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/R1/Axiom_MNEc2M_{}_Annotation.r1.csv".format(design),
            sep=',',
            skiprows=skipno("/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/CallTwo/R1/Axiom_MNEc2M_{}_Annotation.r1.csv".format(design))
        ) 
        metrics = metrics.merge(snp_info,how='left',left_on='probeset_id',right_on='Probe Set ID') 
        allMetrics = pd.concat([allMetrics,metrics])
    allMetrics.to_csv("THE_SNP_LIST.tsv")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
