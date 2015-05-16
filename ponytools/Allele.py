class Allele(object):
    def __init__(self):
        pass    
    @staticmethod
    def vcf2geno(allele_code):
        return {
            './.' : -1,
            '0/0' : 0,
            '0/1' : 1,
            '1/0' : 1,
            '1/1' : 2,
            '2/2' : -1,
            '0/2' : -1,
            '2/0' : -1,
            '1/2' : -1,
            '2/1' : -1
        }[allele_code.replace('|','/')]
    @staticmethod
    def geno2vcf(allele_code):
        return {
            -1 :'./.',
            0  : '0/0',
            1  :'0/1',
            1  :'1/0',
            2  :'1/1',
        }[allele_code]
 
