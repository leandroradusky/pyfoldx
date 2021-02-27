'''
Created on Mar 26, 2013

@author: leandro
'''
from numpy.lib.scimath import sqrt

class Atom(object):
    '''
    classdocs
    '''

    def __init__(self,code,pos,x,y,z,residue = None, bfactor ="1.00"):
        '''
        Constructor
        '''
        self.code = code
        self.pos = pos
        self.x = x
        self.y = y
        self.z = z
        self.residue = residue
        self.bfactor = bfactor
        
    def getInPDBFormat(self, asHetatm=False):
        atom_type = "ATOM  " if not asHetatm else "HETATM"
        atom_pos = "{:5}".format(int(self.pos))
        atom_code = ("  "+self.code+"     ")[0:5]
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
        
    # @param other: another atom
    # @return: distance between atoms 
    def __sub__(self, other):
        return sqrt(pow(self.x - other.x, 2) + pow(self.y - other.y, 2) + pow(self.z - other.z, 2))

    def __eq__(self,other):
        return  self.code == other.code and self.residue==other.residue

    def __repr__(self):
        return "Atom("+self.code+","+str(self.pos)+","+str(self.x)+","+str(self.y)+","+str(self.z)+","+str(self.residue)+","+str(self.bfactor)+")"

    def __str__(self):
        return "Atom("+self.code+","+str(self.pos)+","+str(self.x)+","+str(self.y)+","+str(self.z)+","+str(self.residue)+","+str(self.bfactor)+")"
