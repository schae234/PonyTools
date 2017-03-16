class Chromosome(object): 
    '''
    A Chromosome is a lightweight object which maps indices to 
    string positions.
    '''
    def __init__(self,seq): 
        self.seq = str(seq) 
    def __getitem__(self,pos): 
        if isinstance(pos,slice):
            return self.seq[pos.start-1:pos.stop]
        #chromosomes start at 1, python strings start at 0 
        return self.seq[int(pos)-1]

    def __len__(self):
        return len(self.seq)
 
