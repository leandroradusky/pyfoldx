'''
Created on Nov 2, 2020

@author: lradusky
'''
from pyfoldx.structure.CifParser import parseCifFile
from pyfoldx.structure.PdbParser import parsePdbFile
from pyfoldx.structure.Atom import Atom
from pyfoldx.structure.Residue import Residue
from pyfoldx.handlers.URLRetrieveHandler import URLRetrieveHandler
from pyfoldx.handlers.FileHandler import FileHandler
from pyfoldx.foldx import foldx
from pyfoldx.structure import Ensemble
from pyfoldx.structure.misc import PDB_PATH, ThreeOne

from _collections import defaultdict
import pandas as pd
import re

class Structure(object):
    '''
    Class representing structures (single model)
    
    :ivar data:  The atomic data, dict[chain] -> list of Residue objects 
    '''
    
    def __init__(self, code, path="", from_string=""):
        '''
        Constructor of a structure, will be instantiated from code, file or string depending on its parameters.
        
        :param code: Code of the structure, if no path or string specified this code will be used to download the structure from PDB.
        :param path: (Optional) path of the file to create the structure from, in CIF or PDB format.
        :param from_string: (Optional) String or list of strings as lines to instantiate structure from.
        '''
       
        self.data = defaultdict(lambda:defaultdict(lambda:None))
        self.code = code
        
        if path == "" and from_string == "":
            path = PDB_PATH+code[1]+code[2]+"/"+code+".pdb"
            lines = URLRetrieveHandler.RetrieveFileLines('http://files.rcsb.org/download/'+code+'.pdb')
            FileHandler.writeLines(path, lines)
            parsePdbFile(lines, self)
        if from_string != "":
            parsePdbFile(from_string, self)
        elif path.find(".cif") != -1 or path.find(".cif.gz") != -1:
            parseCifFile(code, path, self)
        elif path.find(".pdb") != -1 or path.find(".pdb.gz") != -1 or \
             path.find(".ent") != -1 or path.find(".ent.gz") != -1:
            lines = FileHandler.getLines(path)
            parsePdbFile(lines, self)
        
        
    def _addPdbAtom(self,fields):
        '''
        Appends an atom to the structure data dictionary based on its properties-
        
        :param fields: List of strings in the same order than a PDB file.
        '''
        [atom_num, atom_code,chain,seqno,resname,x,y,z,bfactor,hetatm] = fields
        
        seqno =  int("".join([s for s in seqno if s.isdigit()]))
        
        if self.data[chain][int(seqno)] is None:
            self.data[chain][int(seqno)] = Residue(resname, int(seqno),chain, [],hetatm)
        
        res = self.data[chain][int(seqno)]
        at = Atom(atom_code,atom_num, float(x),float(y),float(z), res, float(bfactor))
        res.addAtom(at)
    
    def _addCifAtom(self,fields):
        '''
        Appends an atom to the structure data dictionary based on its properties-
        
        :param fields: List of strings in the same order than a CIF file.
        '''
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
        '''
        Appends an atom to the structure data dictionary based on its properties-
        
        :param fields: List of strings in the same order than a file in the format specified in format parameter.
        :param format: Format of the fields passed as list of strings.
        '''
        if format == "PDB":
            self._addPdbAtom(fields)
        elif format == "CIF":
            self._addCifAtom(fields)
    
    def toPdb(self, chain=""):
        '''
        Structure in PDB format.
        
        :param chain: (optional) The chain to return the lines from.
        :return: The lines corresponding to the structure in pdb format.
        '''
        ret = []
        for ch in self.data:
            for res in sorted(self.data[ch].keys()):
                if chain != "" and self.data[ch][res].chain != chain: continue
                for atom in self.data[ch][res].atoms:
                    ret.append(atom.toPdb(self.data[ch][res].hetatm))
        return ret
    
    def getSequence(self, chain):
        '''
        Sequence of a chain within a Structure object.
        
        :param chain: The chain to obtain the sequence from.
        :return: The lines corresponding to the structure in pdb format.
        '''
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
    def getTotalEnergy(self, terms=[],consider_waters=False):
        '''
        Computes FoldX total energy for a structure.
        
        :param terms: List of energy terms to return, if not specified, all the energy terms will be returned.
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing the energy terms computed.
        '''
        print( "Computing total energy for structure..." )
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=foldx.ENERGY_TERMS)
            terms=foldx.ENERGY_TERMS
        
        retDf.loc[self.code] = \
              foldx.getTotalEnergy(self.code, self.toPdb(), consider_waters)[terms].loc[self.code]
              
        print( "Energy computed." )
        return retDf
    
    def getResiduesEnergy(self, considerWaters=False):
        '''
        Computes FoldX energy for all the residues within a Structure object.
        
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing energy terms as columns and residues as index.
        '''
        
        print( "Computing residue energy for structure..." )
        
        retDf = foldx.getResiduesEnergy( self.toPdb(), considerWaters)
        
        firstDigit = lambda s: [int(i) for i in range(0, len(s)) if s[i].isdigit()][0]
        firstPart = lambda x: x[firstDigit(x):]
        secondPart = lambda x: x[0:firstDigit(x)]
        
        retDf.index = retDf.index.map(lambda x: (("000000"+firstPart(x))[-5:]) +"_"+ secondPart(x))
        retDf.sort_index(inplace=True)
        
        print( "Energy computed." )
        return retDf
        
    def repair(self, fix_residues=[], verbose=True):
        '''
        Minimize and complete sidechains of a Structure with FoldX.
        
        :param fix_residues: List of residues in FoldX format to be fixed during repair.
        :return: Structure object with the repaired structure.
        '''
        
        if verbose:
            print( "Repairing structure..." )
        
        ret = Structure(self.code+"_Rep", from_string=foldx.repair(self.toPdb(), fix_residues))
        
        if verbose:
            print( "Structure repaired." )
        return ret
        
    def mutate(self, mutations, number_of_runs=1, terms=[], generate_mutations_ensemble=False, verbose=True):
        '''
        Generate mutated models with FoldX.
        
        :param mutations: List of mutations in FoldX format to be modeled.
        :param number_of_runs: Number of models to be generated.
        :param terms: List of energy terms to return, if not specified, all the energy terms will be returned.
        :param generate_mutations_ensemble: If True, an Ensemble object with the mutations and its wild types is returned.
        :return: Tuple with DataFrame of the ddG of the generated models and ensemble of the mutations if this parameter its true.
        '''
        if verbose:
            print( "Computing mutation(s) %s on target structure..." % mutations )
        
        # if generate mutated ensemble, generate empty ensembles to load
        if generate_mutations_ensemble:
            trajMut = Ensemble.Ensemble(self.code+"_Mut",working_path="")
            trajWT = Ensemble.Ensemble(self.code+"_WT",working_path="")
        else:
            trajMut = None
            trajWT = None
        
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=foldx.ENERGY_TERMS)
            terms=foldx.ENERGY_TERMS

        ddGsDf, mutModels = foldx.mutate(self.toPdb(), mutations, number_of_runs)
        
        # if generate mutated ensemble, add models to the ensemples
        if generate_mutations_ensemble:
            for i in range(number_of_runs):
                trajMut.addFrame(self.code+"_"+ mutations.replace(",","_").replace(";","_"),
                                 "", Structure(self.code+"_"+ mutations.replace(",","_").replace(";","_"), from_string=mutModels[i][0].split("\n")))
                trajWT.addFrame(self.code+"_"+ mutations.replace(",","_").replace(";","_"),
                                 "", Structure(self.code+"_"+ mutations.replace(",","_").replace(";","_"), from_string=mutModels[i][1].split("\n")))
                
                retDf.loc[self.code+"_"+str(i)] = ddGsDf.loc[i]
        
        if verbose:
            print( "Energy computed." )
        return retDf, trajMut, trajWT

if __name__ == "__main__":
    
    print ("started")
    
    st=Structure("2CLD")
    
    print( st.getTotalEnergy() )
    print( st.getResiduesEnergy() )
    
    for line in st.toPdb():
        print( line )
    
    print ("done")
    

