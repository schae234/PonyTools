class Trio(object):
    genomap = { # shared genomap
        './.' : -1, '0/0' : 0, '0/1' : 1, '1/0' : 1, '1/1' : 2,
        -1 : -1, 0 : 0, 1 : 1, 2 : 2     
    }
    def __init__(self,offspring=None,father=None,mother=None):
        self.off = self.genomap[offspring]
        self.fat = self.genomap[father]
        self.mot = self.genomap[mother]

    def consistent(self):
        if self.off == -1: # give benefit of doubt
            consistent = True
        elif self.off == 0 and (self.mot == 0 or self.mot ==1) and (self.fat == 0 or self.fat == 1):
            # both self.father and self.mother must have 0 allele
            consistent = True
        elif self.off == 1 and (not (self.mot == 0 and self.fat == 0) or not (self.mot == 2 and self.fat == 2)):
            # both self.father and self.mother cannot be homozygous
            consistent = True
        elif self.off == 2 and (self.mot == 1 or self.mot ==2) and (self.fat == 1 or self.fat == 2): 
            # both self.father and self.mother must have 1 allele
            consistent = True
        elif self.off == 0 and ((self.mot == -1 and self.fat != 2) or (self.mot != 2 and self.fat == -1)):
            # if one parent is missing, the other must not be opposite homozygous
            consistent = True
        elif self.off == 2 and ((self.mot == -1 and self.fat != 0) or (self.mot != 0 and self.fat == -1)):
            # if one parent is missing, the other must not be opposite homozygous
            consistent = True
        elif self.off == 1 and (self.mot == -1 or self.fat == -1):
            # when self.offspring is homozygous, all bets are self.off
            consistent = True
        else:
            consistent = False
        return consistent
