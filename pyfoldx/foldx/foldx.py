'''
Created on Nov 9, 2020
@author: lradusky
@summary: python bindings for foldx commands
'''

from pyfoldx.structure.misc import OutputGrabber, in_notebook
from pyfoldx.foldx import foldxHandler
from pyfoldx.structure import structure
import pandas as pd
from pyfoldx.foldx.foldxHandler import ENERGY_TERMS
    


def getInterfaceEnergy(pdb_string, consider_waters=False):
    """
    Compute the interface energy of FoldX for molecules within a structure
    
    :param pdb_string: a string containing a structure in PDB format
    :param consider_waters: take waters into account for energy computations
    
    :return: pandas Dataframe with energy terms of the analyzed structure
             with index corresponding to pair of molecules
    """
    
    if type(pdb_string) == type([]):
        pdb_string = "\n".join(pdb_string)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            df = pd.DataFrame.from_dict(foldxHandler.getComplexEnergy(pdb_string, consider_waters), orient="index", columns=ENERGY_TERMS)
    else:
            df = pd.DataFrame.from_dict(foldxHandler.getComplexEnergy(pdb_string, consider_waters), orient="index", columns=ENERGY_TERMS)
    
    return df

def getResiduesEnergy(pdb_string, consider_waters=False):
    """
    Compute the interface energy of FoldX for molecules within a structure
    
    :param pdb_string: a string containing a structure in PDB format
    :param consider_waters: take waters into account for energy computations
    
    :return: pandas Dataframe with energy terms of the analyzed structure
             with index corresponding to each residue
    """
    
    if type(pdb_string) == type([]):
        pdb_string = "\n".join(pdb_string)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            df = pd.DataFrame.from_dict(foldxHandler.getResiduesEnergy(pdb_string, consider_waters), orient="index", columns=ENERGY_TERMS)
    else:
        df = pd.DataFrame.from_dict(foldxHandler.getResiduesEnergy(pdb_string, consider_waters), orient="index", columns=ENERGY_TERMS)
    
    return df

def mutate(pdb_string, mutations, number_of_runs=1):
    """
    Generate a mutated structure
    
    :param pdb_string: a string containing a structure in PDB format
    :param mutations: mutations to be performed in foldx format (example GI1A means to mutate GLY in position 1 of molecule I to ALA)
    
    :return: tuple, where element 0 is the DataFrame of ddGs of each models respect to wildtype structure 
            and element 1 is the array of strings with the mutated model
    """
    
    if type(pdb_string) == type([]):
        pdb_string = "\n".join(pdb_string)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            ddGs, mutModels = foldxHandler.getMutant(pdb_string, mutations, number_of_runs) 
    else:
        ddGs, mutModels = foldxHandler.getMutant(pdb_string, mutations, number_of_runs) 
    
    return (pd.DataFrame(ddGs, columns=ENERGY_TERMS), mutModels)
    
def repair(pdb_string,fix_residues=[]):
    """
    Repair the sidechains of a structure
    
    :param pdb_string: a string containing a structure in PDB format
    :param fix_residues: list of residues to remain fixed in foldx format (example GI1 means GLY in position 1 of molecule I)
    
    :return: repaired model in string format
    """
    
    if type(pdb_string) == type([]):
        pdb_string = "\n".join(pdb_string)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            ret = foldxHandler.getRepairedStructure(pdb_string, fix_residues)
    else:
        ret = foldxHandler.getRepairedStructure(pdb_string, fix_residues)
    
    return ret

def getNetworks(pdb_string):
    
    """
    Get networking information of a structure
    
    :param pdb_string: a string containing a structure in PDB format

    :return: dictionary where keys are network types and elements are linked residue names for that network
    """
    
    if type(pdb_string) == type([]):
        pdb_string = "\n".join(pdb_string)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            ret = foldxHandler.getNetworks(pdb_string)
    else:
        ret = foldxHandler.getNetworks(pdb_string)
    
    return ret

if __name__ == "__main__":
    print("started")

    code = "6U6M"
    st = structure(code)
    print( getTotalEnergy(code, st.toPdb(), False) )
    print( getInterfaceEnergy(st.toPdb(), False) )
    #print( getResiduesEnergy(st.toPdb(), False) )
    #print( mutate(st.toPdb(), "GI29A;",1)[0] )
    #print( mutate(st.toPdb(), "GI29A;",1)[1][0][0] )
    
    #print( repair(st.toPdb(), ["GI29"]) )
    #for k,v in getNetworks(st.toPdb()).items() : print( k , v )
    '''
    code = "5XJL"
    st = structure(code)
    print (getInterfaceEnergy(st.toPdb(), False) )
    '''
    
    
    
    
    