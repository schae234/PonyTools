import pandas as pd
import os, errno
import sys
from subprocess import Popen, PIPE
from optparse import OptionParser

def skipno(filename):
    with open(filename,'r') as IN:
        skip = 0
        while IN.readline().startswith("#"):
            skip += 1
        return skip


def main(args):

    parser = OptionParser()
    options,files = parser.parse_args() 

    # Read in samples
    samples = pd.read_csv("2MequineAX_Sample_Table.csv")
    # Read in sample type
    sample_type = pd.read_table("SampleType.csv")
    samples = samples.merge(sample_type[['SampleName','Isolatedfrom']].drop_duplicates(),how="left",left_on="SampleName",right_on="SampleName")
    
    
    # no samples should have a NULL method oisolation
    assert(any(samples.Isolatedfrom.isnull()) == False)
    
    # Generate different CEL lists for different sampe types:
    cel_path = "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/raw_data"
    
    processes = []
    for groups,df in samples.groupby(["Design","Isolatedfrom"]):
        design,tissuetype = groups
        path = os.path.join(*groups)
        # make the batch directory
        os.makedirs(path,exist_ok = True)
        # Make the CELs file
        with open(os.path.join(path,"CELS.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            # print out the cels belonging to that batch
            for CEL in df.BestArray.values:
                print(os.path.join(cel_path,str(CEL)+".CEL"),file=CELS)
        # EXECUTE!
        stdout = open(os.path.join(path,'stdout.txt'),'w')
        stderr = open(os.path.join(path,'stderr.txt'),'w')
        # Dont run script if the output is already there
        if not os.path.exists(os.path.join(path,"qc-geno.txt")):
            processes.append(Popen([
                'apt-geno-qc',
                '--analysis-files-path', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/'.format(design,design),
                '--xml-file', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/Axiom_MNEc2M_{}.v1.apt-geno-qc.AxiomQC1.xml'.format(design,design,design),
                '--cel-files', os.path.join(path,'CELS.txt'),
                '--out-dir', os.path.join(path), 
                '--out-file', os.path.join(path,'qc-geno.txt'),
                '--log-file', os.path.join(path,'log.txt')
            ],stdout=stdout,stderr=stderr))
    # Wait for apt-geno-qc to finish
    while processes:
        for p in processes[:]:
            if p.poll() is not None:
                print("Process {} ended with exit code {}".format(p.pid,p.returncode))
                processes.remove(p)
    
    # Birdseed needs to know about SNPs with fucked up shit
    special_snps = pd.DataFrame()
    special_path = '/project/mccuelab/rob/MNEc2M/special_snps.txt'
    for i,design in enumerate(['A','B','C']):
        annot = "/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Annotation.v1/CsvAnnotationFile_{}.v1.txt".format(design,design,i+1)
        special_snps = pd.concat([special_snps,pd.read_csv(annot,skiprows=skipno(annot),sep=',',quotechar='"')])
    special_snps = special_snps[[x in ['chrX','chrY'] for x in special_snps.cust_chr]]
    special_snps = special_snps[["Affy SNP ID","cust_chr"]]
    special_snps.columns = ['probeset_id','chr']
    special_snps['copy_female'] = [2 if chr == 'chrX' else 0  for chr in special_snps.chr ]
    special_snps['copy_male'] = 1
    special_snps.to_csv(special_path,sep="\t",index=False)
    
    ALLDQC = pd.DataFrame()
    # Remove Samples with DQC less than 0.6
    for groups,df in samples.groupby(['Design','Isolatedfrom']):
        design,tissuetype = groups
        path = os.path.join(*groups)
        # Read in the qc-geno results and filter out DQC below 0.6
        DQCResults = pd.read_table(os.path.join(path,"qc-geno.txt"),skiprows=skipno(os.path.join(path,"qc-geno.txt")))
        DQCResults['Design'] = design
        DQCResults['TissueType'] = tissuetype
        ALLDQC = pd.concat([ALLDQC,DQCResults])
        with open(os.path.join(path,"CELS2.txt"),'w') as CELS:
            print("cel_files",file=CELS)
            for i,CEL in DQCResults.iterrows():
                if CEL.axiom_dishqc_DQC < 0.6:
                    print("Dropping {} for {},{} because it was < 0.6".format(CEL.cel_files,design,tissuetype))
                else:
                    print(os.path.join(cel_path,CEL.cel_files),file=CELS)
        # Produce Genotypes
        if not os.path.exists(os.path.join(path,"AxiomGT1.calls.txt")):
            os.makedirs(os.path.join(path,'calls'),exist_ok = True)
            stdout = open(os.path.join(path,'calls','geno-stdout.txt'),'w')
            stderr = open(os.path.join(path,'calls','geno-stderr.txt'),'w')
            processes.append(Popen([
                'apt-probeset-genotype',
                '--analysis', 'brlmm-p',
            #'--write-models',
            #"--no-gender-force",
                '--analysis-files-path', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/'.format(design,design),
                "--cdf-file"           ,"Axiom_MNEc2M_{}.v1.cdf".format(design),
            #'--analysis', 'birdseed-dev',
            #'--analysis-files-path', '/project/mccuelab/rob/MNEc2M/Affy_MNEc2M/cluster_{}/AffyCalls/Axiom_MNEc2M_{}_Analysis.v1/'.format(design,design),
                '--set-analysis-name','AxiomGT1-BS',
                '--chip-type', 'Axiom_MNEc2M_{}'.format(design),
                '--read-models-brlmmp', os.path.join(path,'brlmm-p-plus.force.snp-posteriors.txt'),
            #'--summaries',
                "--chrX-probes"        ,"Axiom_MNEc2M_{}.v1.chrXprobes".format(design),
            #"--special-snps"       ,"/project/mccuelab/rob/MNEc2M/special_snps.txt",
                "--chrY-probes"        ,"Axiom_MNEc2M_{}.v1.chrYprobes".format(design),
            #"--target-sketch"      ,"Axiom_MNEc2M_{}.v1.AxiomGT1.sketch".format(design),
                "--set-gender-method"  ,"cn-probe-chrXY-ratio",
                "--em-gender"          ,"false",
                "--female-thresh"      ,"0.5",
                "--male-thresh"        ,"0.9", 
                '--cel-files', os.path.join(path,'CELS2.txt'),
                '--out-dir', os.path.join(path,"calls"), 
                '--log-file', os.path.join(path,'calls','genotype-log.txt')
            ],stdout=stdout,stderr=stderr))
    while processes:
        for p in processes[:]:
            if p.poll() is not None:
                print("Process {} ended with exit code {}".format(p.pid,p.returncode))
                processes.remove(p)
    
    

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
