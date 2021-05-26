'''
Created on Nov 2, 2020

@author: lradusky
'''
from pyfoldx.structure.cifParser import parseCifFile
from pyfoldx.structure.pdbParser import parsePdbFile
from pyfoldx.structure.atom import Atom
from pyfoldx.structure.residue import Residue
from pyfoldx.handlers.urlRetrieveHandler import URLRetrieveHandler
from pyfoldx.handlers.fileHandler import FileHandler
from pyfoldx.foldx import foldxHandler
from pyfoldx.structure.ensemble import Ensemble
from pyfoldx.structure.misc import PDB_PATH, ThreeOne

from _collections import defaultdict
import pandas as pd
import re

class Structure(object):
    '''
    Class representing structures (single model)
    
    :ivar code: The code given to the structure (for loading and saving) 
    :ivar data: The atomic data, dict[chain] -> list of residue objects 
    :ivar resolution: If parsed from full file, structure resolution is parsed and set. 
    
    :ivar totalEnergy: pandas DataFrame with total energy terms
    :ivar interfaceEnergy: pandas DataFrame with total energy terms
    :ivar residuesEnergy: pandas DataFrame with energy terms per residue
    :ivar molEnergy: pandas DataFrame with energy per molecule computed when complex energy is called
    
    :ivar networks: dictionary of DataFrames with networks computed
    :ivar interfaceResidues: dict of pair of molecules and their interface residues
    :ivar unrecognizedMolecules: residues not recognized by FoldX
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
        self.resolution = None
        
        # FoldX attributes, setted upon calling to getters
        self.unrecognizedMolecules = []
        self.interfaceResidues = None
        
        self.totalEnergy = None
        self.interfaceEnergy = None
        self.molEnergy = None
        self.residuesEnergy = None
        self.networks = None
        
        if path == "" and from_string == "":
            path = PDB_PATH+code[1]+code[2]+"/"+code+".pdb"
            lines = URLRetrieveHandler.RetrieveFileLines('http://files.rcsb.org/download/'+code+'.pdb')
            FileHandler.writeLines(path, lines)
            parsePdbFile(lines, self)
            self.resolution = self.setResolution(lines)
        if from_string != "":
            parsePdbFile(from_string, self)
        elif path.find(".cif") != -1 or path.find(".cif.gz") != -1:
            parseCifFile(code, path, self)
        elif path.find(".pdb") != -1 or path.find(".pdb.gz") != -1 or \
             path.find(".ent") != -1 or path.find(".ent.gz") != -1:
            lines = FileHandler.getLines(path)
            parsePdbFile(lines, self)
            self.resolution = self.setResolution(lines)
        
        
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
    
    def toPdbFile(self, path, chain=""):
        '''
        Save structure in PDB format.
        
        :param path: path of the file to save to.
        :param chain: (optional) The chain to return the lines from.
        '''
        FileHandler.writeLines( path, self.toPdb(chain) )
    
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
    
    def getMinMaxBFactor(self):
        '''
        Minimum and maximum bfactor to normalize residue bfactors
        
        :return: Min and Max bfactors
        '''
        mi = -1
        ma = 10000
        for ch in self.data:
            for res in sorted(self.data[ch].keys()):
                bf = self.data[ch][res].getAvgBFactor()
                if bf < ma:
                    ma = bf
                if bf > mi:
                    mi = bf
        
        return ma, mi
    
    def setResolution(self, lines):
        '''
        Return resolution of the crystal if defined
        
        :return: Resolution of the crystal
        '''
        try:
            for line in lines:
                if line.find('REMARK   2 RESOLUTION.') != -1:
                    return float(line[22:30].strip())
        except Exception as e:
            return None

    ##################
    # FoldX Operations
    ##################
    def getTotalEnergy(self, terms=[],consider_waters=False):
        '''
        Computes FoldX total energy for a structure.
        
        :param terms: List of energy terms to return, if not specified, all the energy terms will be returned.
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing the energy terms computed. Total energy attribute is setted to this dataframe as well.
        '''
        print( "Computing total energy for structure..." )
        if terms != []:
            self.totalEnergy = pd.DataFrame(columns=terms)
        else:
            self.totalEnergy = pd.DataFrame(columns=foldxHandler.ENERGY_TERMS)
            terms=foldxHandler.ENERGY_TERMS
        
        foldxHandler.getTotalEnergy(self, consider_waters=consider_waters)[terms].loc[self.code]

        
        print( "Energy computed." )
        return self.totalEnergy
    
    def getInterfaceEnergy(self,consider_waters=False, verbose=True):
        '''
        Computes FoldX total energy for a structure.
        
        :param terms: List of energy terms to return, if not specified, all the energy terms will be returned.
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing the energy terms computed by pairs of molecules. Also setted in interfaceEnergy attribute. Other attributes setted:
                    individual_energies: DataFrame, of each molecule.
                    interface_residues: List, between pairs of molecules.
                    interaction_energies: Detailed energies term by term, DataFrame
                    summary: Returned dataframe
                    unrecognized_molecules: List
        '''
        if verbose:
            print( "Computing complex energy for structure..." )
        
        foldxHandler.getComplexEnergy( self, consider_waters )
              
        if verbose:
            print( "Energy computed." )
        return self.interfaceEnergy
    
    def getResiduesEnergy(self, consider_waters=False):
        '''
        Computes FoldX energy for all the residues within a Structure object.
        
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing energy terms as columns and residues as index.
        '''
        
        print( "Computing residue energy for structure..." )
        
        foldxHandler.getResiduesEnergy( self, consider_waters)
        
        print( "Energy computed." )
        return self.residuesEnergy
        
    def getNetworks(self, verbose=True):
        '''
        Computes FoldX energy for all the residues within a Structure object.
        
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing energy terms as columns and residues as index.
        '''
        if verbose:
            print( "Computing networks within structure..." )
        
        foldxHandler.getNetworks( self )
        
        if verbose:
            print( "Networks computed." )
    
    def alanineScan(self, consider_waters = False, verbose=True):
        '''
        Computes FoldX energy upon mutation to alanine for all the residues within the structure
        
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing the energy variation upon mutation to alanine for all the residues within the structure.
        '''
        if verbose:
            print( "Performing Alanine Scan..." )
            
        ret = foldxHandler.getAlaScan(self,consider_waters)
        
        if verbose:
            print( "Alanine Scan finished." )
            
        return ret
    
    def positionScan(self, positions, chain, consider_waters = False, verbose=True):
        '''
        Computes FoldX  position specific scoring matrix (PSSM) for specified residues within the structure
        
        :param positions: list of positions to be scanned
        :param chain: chain to be considered for specified positions
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing the energy variation upon mutation to alanine for all the residues within the structure.
        '''
        if verbose:
            print( "Performing PSSM..." )
        
        positions_string = ""
        for pos in positions:
            rescode = ThreeOne[self.data[chain][pos].code]
            positions_string += rescode+chain+str(pos)+"a,"
        positions_string = positions_string[0:-1]
        
        ret = foldxHandler.getPSSM(self, positions_string)
        
        if verbose:
            print( "PSSM finished." )
            
        return ret
    
    def repair(self, fix_residues=[], verbose=True):
        '''
        Minimize and complete sidechains of a Structure with FoldX.
        
        :param fix_residues: List of residues in FoldX format to be fixed during repair.
        :return: Structure object with the repaired structure.
        '''
        
        if verbose:
            print( "Repairing structure..." )
        
        ret = foldxHandler.getRepairedStructure(self, fix_residues)
        
        if verbose:
            print( "Structure repaired." )
        return ret
        
    def mutate(self, mutations, number_of_runs=1, terms=[], verbose=True):
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
        trajMut = Ensemble(self.code+"_Mut",working_path="")
        trajWT = Ensemble(self.code+"_WT",working_path="")
        
        ddGsDf, mutModels, wtModels = foldxHandler.getMutants(self, mutations, number_of_runs)
        
        # if generate mutated ensemble, add models to the ensemples
        for i in range(number_of_runs):
            trajMut.addFrame(mutModels[i].code, "", mutModels[i])
            trajWT.addFrame(wtModels[i].code, "", wtModels[i])
    
        if verbose:
            print( "Energy computed." )
            
        return ddGsDf, trajMut, trajWT

if __name__ == "__main__":
    
    print ("started")
    
    code = "2CI2"
    st = Structure(code)
    print( st.alanineScan() )
    print( st.positionScan(range(75,77),"I") )
    
    print( st.getTotalEnergy() )
    repSt = st.repair()
    print( repSt.getTotalEnergy() )
    st.getNetworks()
    print( st.networks["Electro"] )
      
    code = "1KX1"
    st = Structure(code)
    print( st.getInterfaceEnergy() )
    print( st.interfaceResidues )
    print( st.molEnergy )
     
    print( st.getResiduesEnergy() )
    ddGsDf, trajMut, trajWT = st.mutate( "AA73G;",3 )
    print( ddGsDf )
    print( trajMut )
    print( trajMut.getFrame(0).totalEnergy )
    
    st = Structure("2B3B")
    st_A = Structure( "2B3B_A", from_string= st.toPdb(chain="A") )
    energyDf_GLC, ensMut_GLC, ensWT_GLC          = st_A.mutate("DA278H;")
    
     
    print ("done")
    

