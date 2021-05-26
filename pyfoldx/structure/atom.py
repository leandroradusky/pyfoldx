'''
Created on Mar 26, 2013

:author: lradusky
'''

from numpy import sqrt

class Atom(object):
    '''
    Class representing atoms within structures.
    '''

    def __init__(self,code,pos,x,y,z,residue = None, bfactor ="1.00"):
        '''
        Constructor of Atom, with all the atributes they have in a structure
        '''
        
        self.code = code
        self.pos = pos
        self.x = x
        self.y = y
        self.z = z
        self.residue = residue
        self.bfactor = bfactor
        
    def toPdb(self, as_hetatm=False):
        '''
        A string for the atom object in pdb format.
        
        :param as_hetatm: The returned line is identified as heteroatom.
        :return: The line corresponding to the atom in pdb format.
        '''
        
        atom_type = "ATOM  " if not as_hetatm else "HETATM"
        atom_pos = "{:5}".format(int(self.pos))
        atom_code = ("  "+self.code.replace("'","*")+"     ")[0:5]
        atom_res_code = ("      "+self.residue.code)[-4:]
        atom_res_chain = ("      "+self.residue.chain)[-2:]
        atom_res_pos = ("      "+str(self.residue.pos))[-4:]
        atom_x = "{:12.3f}".format(self.x)
        atom_y = "{:8.3f}".format(self.y)
        atom_z = "{:8.3f}".format(self.z)
        atom_ocuppancy = "  1.00"
        atom_bfactor = "{:6.2f}".format(float(self.bfactor))
        atom_tail = "           "+self.code[0]+"  "
        
        return atom_type+atom_pos+atom_code+atom_res_code+atom_res_chain+atom_res_pos+\
                atom_x+atom_y+atom_z+atom_ocuppancy+atom_bfactor+atom_tail
        
    def __sub__(self, other):
        '''
        Substraction of two atom vectors, returning its distance.
        
        :param other: another atom to measure distance with.
        :return: distance between self atom and other athom. 
        '''
        return sqrt(pow(self.x - other.x, 2) + pow(self.y - other.y, 2) + pow(self.z - other.z, 2))

    def __eq__(self,other):
        '''
        Equal operator check if two atoms have the same atom code within the same Residue.
        
        :param other: another atom to measure compare with.
        :return: True if they are the same atom code within the same Residue. 
        '''
        return  self.code == other.code and self.residue==other.residue

    def __repr__(self):
        '''
        Representation of the atom.
        
        :return: String with the atom atributes.
        '''
        return self.__str__()

    def __str__(self):
        '''
        Representation of the atom.
        
        :return: String with the atom atributes.
        '''
        return "Atom("+self.code+","+str(self.pos)+","+str(self.x)+","+str(self.y)+","+str(self.z)+","+str(self.residue)+","+str(self.bfactor)+")"




