'''
Created on Nov 2, 2020
@author: lradusky
'''
from src.handlers.FileHandler import FileHandler

def parsePdbFile(lines, struct):
    
    if (type(lines)==type("")):
        lines = lines.split("\n")
    
    for line in lines:
        if (line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM' ):
            hetatm=False
            if line[0:6].strip() == 'HETATM':
                hetatm=True
            # Extract the important fields
            atom_num = line[6:11].strip()
            atom_code = line[11:16].strip()
            chain = line[21:22].strip()
            seqno =line[22:27].strip()
            resname = line[17:20].strip()
            try: x = float(line[30:38].strip())
            except: x = 0.0
            try: y = float(line[38:46].strip())
            except: y = 0.0
            try: z = float(line[46:54].strip())
            except: z = 0.0
            try: bfactor = float(line[60:66].strip())
            except: bfactor = 0.0
            
            struct.addAtom([atom_num, atom_code,chain,seqno,resname,x,y,z,bfactor,hetatm], "PDB")