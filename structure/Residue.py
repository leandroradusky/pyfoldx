'''
Created on Mar 25, 2013

@author: lradusky
'''

class Residue(object):
    '''
    Class representing residues within structures
    '''

    # @param code: three letter code of the residue
    # @param pos: position of the residue in the chain
    # @param atoms: list of atoms of the residue   
    # @param hetatm: is hetatm 
    def __init__(self, code, pos, chain = None, atoms = [], hetatm = False):
        '''
        Constructor of Residue, with all the atributes they have in a structure and a list of its atoms.
        
        :param code: The residue code (three letters)
        :param pos: The residue number.
        :param chain: The chain the residue belongs to.
        :param atoms: List of atom objects belonging this residue.
        :param hetatm: True if the residue is an ligand, False if its a residue. 
        '''
        self.code = code
        self.pos = pos
        self.atoms = atoms
        self.hetatm = hetatm
        self.chain = chain
    
    def addAtom(self, atom):
        '''
        Appends an Atom object to the list of atoms if it is not.
        
        :param atom: Atom object to be included.
        '''
        
        if atom not in self.atoms:
            self.atoms += [atom]
    
    def toPdb(self):
        '''
        String for the atoms of the residue in pdb format.
        
        :return: The lines corresponding to the residue in pdb format.
        '''
        return "\n".join([a.toPdb(self.hetatm) for a in self.atoms])
    
    # @return: average BFactor of the residue
    def getAvgBFactor(self):
        '''
        Average Beta Factor of a residue within all its atom.
        
        :return: Average Beta Factor of a residue within all its atom.
        '''
        atomC = 0.
        bfacS = 0.
        for atom in self.atoms:
            atomC += 1
            bfacS += atom.bfactor
            
        try:
            return bfacS/atomC
        except:
            return 0.0
        
    def __sub__(self, other):
        '''
        Substraction operator between two residues.
        
        :param other: another residue to measure distance with.
        :return: minimun distance between atoms of self Residue and other Residue.
        '''
        min_dist=-1
        for myatom in self.atoms:
            for hisatom in other.atoms:
                if myatom - hisatom < min_dist or min_dist == -1:
                    min_dist = myatom - hisatom
                    
        return min_dist
    
    def __eq__(self,other):
        '''
        Equal operator check if they have the same code and position.
        
        :param other: another Residue compare with.
        :return: True if self and other have the same code and position.
        '''
        return self.code == other.code and self.pos == other.pos
    
    def __repr__(self):
        '''
        Representation of the Residue.
        
        :return: String with the Residue atributes.
        '''
        return "Residue("+self.code+","+str(self.pos)+","+str(self.chain)+","+str(self.hetatm)+")"
    
    def __str__(self):
        '''
        Representation of the Residue.
        
        :return: String with the Residue atributes.
        '''
        return self.__repr__()




