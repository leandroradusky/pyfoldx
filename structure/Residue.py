'''
Created on Mar 25, 2013

@author: leandro
'''

class Residue(object):
    '''
    classdocs
    '''

    # @param code: three letter code of the residue
    # @param pos: position of the residue in the chain
    # @param atoms: list of atoms of the residue   
    # @param hetatm: is hetatm 
    def __init__(self, code, pos, chain = None, atoms = [], hetatm = False):
        '''
        Constructor
        '''
        
        self.code = code
        self.pos = pos
        self.atoms = atoms
        self.hetatm = hetatm
        self.chain = chain
    
    def addAtom(self, atom):
        self.atoms += [atom]
    
    def getInPDBFormat(self):
        return "\n".join([a.getInPDBFormat(self.hetatm) for a in self.atoms])
    
    # @return: average BFactor of the residue
    def getAvgBFactor(self):
        atomC = 0.
        bfacS = 0.
        for atom in self.atoms:
            atomC += 1
            bfacS += atom.bfactor
            
        try:
            return bfacS/atomC
        except:
            return 0.0
        
    # @param other: another residue
    # @return: minimun distance between atoms of residues 
    def __sub__(self, other):
        min_dist=-1
        for myatom in self.atoms:
            for hisatom in other.atoms:
                if myatom - hisatom < min_dist or min_dist == -1:
                    min_dist = myatom - hisatom
                    
        return min_dist
    
    def __eq__(self,other):
        return self.code == other.code and self.pos == other.pos
    
    def __repr__(self):
        return "Residue("+self.code+","+str(self.pos)+","+str(self.chain)+","+str(self.hetatm)+")"
    
    def __str__(self):
        return self.__repr__()




