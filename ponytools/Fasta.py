from collections import defaultdict
from ponytools.Tools import chromosome, log

class Fasta(object):
    def __init__(self):
        self.added_order = []
        self.chroms = {}
        self.attributes = defaultdict(list)

    def add_chrom(self,chrom_name,sequence):
        log("Adding new chromosome: {}",chrom_name)
        self.added_order.append(chrom_name)
        self.chroms[chrom_name] = sequence
    def add_attribute(self,chrom_name,attr):
        self.attributes[chrom_name].append(attr)

    def __getitem__(self,chrom_name):
        return self.chroms[chrom_name]

    @classmethod
    def from_file(cls,fasta_file):
        fasta = cls()
        with open(fasta_file,'r') as IN: 
            cur_chrom = None
            cur_seqs = []
            for line in IN:
                line = line.strip()
                if line.startswith('>'):
                    if cur_chrom:
                        fasta.add_chrom(cur_chrom,chromosome("".join(cur_seqs)))
                    name,*attrs = line.lstrip('>').split()
                    cur_chrom = name
                    cur_seqs = []
                    for attr in attrs:
                        fasta.add_attribute(attr)
                else:
                    cur_seqs.append(line)
            # Add the last chromosome
            fasta.add_chrom(cur_chrom,chromosome("".join(cur_seqs)))
        return fasta
