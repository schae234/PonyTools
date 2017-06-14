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
            '2/1' : -1,
            '1/2' : -1,
            '3/2' : -1,
            '2/3' : -1,
            '3/1' : -1,
            '1/3' : -1
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

    @staticmethod
    def vcf2allele(ref,alt):
        def closure(allele_code):
            return {
                './.' : '--',
                '0/0' : '{}{}'.format(ref,ref),
                '0/1' : '{}{}'.format(alt,ref),
                '1/0' : '{}{}'.format(alt,ref),
                '1/1' : '{}{}'.format(alt,alt),
                '2/2' : -1,
                '0/2' : -1,
                '2/0' : -1,
                '1/2' : -1,
                '2/1' : -1
            }[allele_code.replace('|','/')]
        return closure

