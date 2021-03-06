'''
Created on Nov 2, 2020

@author: lradusky
'''
from src.structure.CifParser import parseCifFile
from src.structure.PdbParser import parsePdbFile
from src.structure.Atom import Atom
from src.structure.Residue import Residue
from src.handlers.URLRetrieveHandler import URLRetrieveHandler
from src.handlers.FileHandler import FileHandler
from src.pyFoldX import pyFoldX
import src.structure.Ensemble as Ensemble

from _collections import defaultdict
import pandas as pd
import re

PDB_PATH = '/data/pdb/divided/pdb/'

ThreeOne = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F'
,'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M'
,'ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T'
,'VAL':'V','TRP':'W','TYR':'Y'}

class Structure(object):

    def __init__(self, code, path="", fromString=""):
        # The atomic data, dict[chain] -> list of residue objects
        self.data = defaultdict(lambda:defaultdict(lambda:None))
        self.code = code
        
        if path == "" and fromString == "":
            path = PDB_PATH+code[1]+code[2]+"/"+code+".pdb"
            lines = URLRetrieveHandler.RetrieveFileLines('http://files.rcsb.org/download/'+code+'.pdb')
            FileHandler.writeLines(path, lines)
            parsePdbFile(lines, self)
        if fromString != "":
            parsePdbFile(fromString, self)
        elif path.find(".cif") != -1 or path.find(".cif.gz") != -1:
            parseCifFile(code, path, self)
        elif path.find(".pdb") != -1 or path.find(".pdb.gz") != -1 or \
             path.find(".ent") != -1 or path.find(".ent.gz") != -1:
            lines = FileHandler.getLines(path)
            parsePdbFile(lines, self)
        
        
    def addPdbAtom(self,fields):
        [atom_num, atom_code,chain,seqno,resname,x,y,z,bfactor,hetatm] = fields
        
        if self.data[chain][int(seqno)] is None:
            self.data[chain][int(seqno)] = Residue(resname, int(seqno),chain, [],hetatm)
        
        res = self.data[chain][int(seqno)]
        at = Atom(atom_code,atom_num, float(x),float(y),float(z), res, float(bfactor))
        res.addAtom(at)
    
    def addCifAtom(self,fields):
        [group_PDB,id,type_symbol,label_atom_id,label_alt_id,label_comp_id,label_asym_id,
        label_entity_id,label_seq_id,pdbx_PDB_ins_code,Cartn_x,Cartn_y,Cartn_z,occupancy,
        B_iso_or_equiv,pdbx_formal_charge,auth_seq_id,auth_comp_id,auth_asym_id,auth_atom_id,
        pdbx_PDB_model_num] = fields
        
        if self.data[auth_asym_id][int(auth_seq_id)] is None:
            self.data[auth_asym_id][int(auth_seq_id)] = Residue(label_comp_id, int(auth_seq_id), \
                                                                 auth_asym_id, [],\
                                                             True if group_PDB=="HETATM" else False)
        
        res = self.data[auth_asym_id][int(auth_seq_id)]
        at = Atom(label_atom_id,id, float(Cartn_x),float(Cartn_y),float(Cartn_z), res, float(B_iso_or_equiv))
        res.addAtom(at)
        
    def addAtom(self, fields, format='PDB'):
        if format == "PDB":
            self.addPdbAtom(fields)
        elif format == "CIF":
            self.addCifAtom(fields)
    
    def toPdb(self, chain=""):
        ret = []
        for ch in self.data:
            for res in sorted(self.data[ch].keys()):
                if chain != "" and self.data[ch][res].chain != chain: continue
                for atom in self.data[ch][res].atoms:
                    ret.append(atom.getInPDBFormat(self.data[ch][res].hetatm))
        return ret
    
    def getSequence(self, chain):
        ret = ""
        for res in self.data[chain]:
            if not self.data[chain][res].hetatm:
                try:
                    ret+=ThreeOne[self.data[chain][res].code]
                except: ret+="-"
        return ret
    
    ##################
    # FoldX Operations
    ##################
    def getTotalEnergy(self, terms=[],considerWaters=False):
        
        print( "Computing total energy for structure..." )
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=pyFoldX.ENERGY_TERMS)
            terms=pyFoldX.ENERGY_TERMS
        
        retDf.loc[self.code] = \
              pyFoldX.getTotalEnergy(self.code, self.toPdb(), considerWaters)[terms].loc[self.code] 

        
        print( "Energy computed." )
        return retDf
    
    def getResiduesEnergy(self, considerWaters=False):
        
        print( "Computing residue energy for structure..." )
        
        retDf = pyFoldX.getResiduesEnergy( self.toPdb(), considerWaters)
        
        firstDigit = lambda s: [int(i) for i in range(0, len(s)) if s[i].isdigit()][0]
        firstPart = lambda x: x[firstDigit(x):]
        secondPart = lambda x: x[0:firstDigit(x)]
        
        retDf.index = retDf.index.map(lambda x: (("000000"+firstPart(x))[-5:]) +"_"+ secondPart(x))
        retDf.sort_index(inplace=True)
        
        print( "Energy computed." )
        return retDf
        
    def repair(self, fixResidues=[], verbose=True):
        
        if verbose:
            print( "Repairing structure..." )
        
        ret = Structure(self.code+"_Rep", fromString=pyFoldX.repair(self.toPdb(), fixResidues))
        
        if verbose:
            print( "Structure repaired." )
        return ret
        
    def mutate(self, mutations, numberOfRuns=1, terms=[], generateMutationsEnsemble=False, verbose=True):
        '''
        :todo: implement several runs and keep the best
        '''
        if verbose:
            print( "Computing mutation(s) %s on target structure..." % mutations )
        
        # if generate mutated ensemble, generate empty ensembles to load
        if generateMutationsEnsemble:
            trajMut = Ensemble.Ensemble(self.code+"_Mut",workingPath="")
            trajWT = Ensemble.Ensemble(self.code+"_WT",workingPath="")
        else:
            trajMut = None
            trajWT = None
        
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=pyFoldX.ENERGY_TERMS)
            terms=pyFoldX.ENERGY_TERMS

        ddGsDf, mutModels = pyFoldX.mutate(self.toPdb(), mutations, numberOfRuns)
        
        # if generate mutated ensemble, add models to the ensemples
        if generateMutationsEnsemble:
            for i in range(numberOfRuns):
                trajMut.addFrame(self.code+"_"+ mutations.replace(",","_").replace(";","_"),
                                 "", Structure(self.code+"_"+ mutations.replace(",","_").replace(";","_"), fromString=mutModels[i][0].split("\n")))
                trajWT.addFrame(self.code+"_"+ mutations.replace(",","_").replace(";","_"),
                                 "", Structure(self.code+"_"+ mutations.replace(",","_").replace(";","_"), fromString=mutModels[i][1].split("\n")))
                
        retDf.loc[self.code+"_"+str(i)] = ddGsDf.loc[i]
        
        if verbose:
            print( "Energy computed." )
        return retDf, trajMut, trajWT

if __name__ == "__main__":
    
    print ("started")
    
    st=Structure("2ci2")
    
    print( st.getTotalEnergy() )
    print( st.getResiduesEnergy() )
    print( st.mutate("GI29A;",numberOfRuns=5,generateMutationsEnsemble=True) )
    
    print ("done")
    

