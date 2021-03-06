'''
Created on Nov 9, 2020
@author: lradusky
@summary: python bindings for foldx commands
'''

import pandas as pd
from src.structure.misc import OutputGrabber, in_notebook
import importlib
import sys
import pyFoldXpp
    
ENERGY_TERMS = ['total','backHbond','sideHbond','energy_VdW','electro','energy_SolvP','energy_SolvH','energy_vdwclash',
                'entrop_sc','entrop_mc','sloop_entropy','mloop_entropy','cis_bond','energy_torsion','backbone_vdwclash',
                'energy_dipole','water','disulfide','energy_kon','partcov','energyIonisation','entr_complex']


def getTotalEnergy(pdbCode, pdbString, considerWaters=False):
    """
    Compute the stability energy of FoldX for a structure
    
    :param pdbCode: pdb code
    :param pdbString: a string containing a structure in PDB format
    :param considerWaters: take waters into account for energy computations
    
    :return: pandas Dataframe with energy terms of the analyzed structure
    """
    
    if type(pdbString) == type([]):
        pdbString = "\n".join(pdbString)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            df= pd.DataFrame(data = [pyFoldXpp.getEnergy(pdbString, considerWaters)],
                         columns= ENERGY_TERMS, index=[pdbCode])
    else:
        df= pd.DataFrame(data = [pyFoldXpp.getEnergy(pdbString, considerWaters)],
                     columns= ENERGY_TERMS, index=[pdbCode])
    
    return df

def getInterfaceEnergy(pdbString, considerWaters=False):
    """
    Compute the interface energy of FoldX for molecules within a structure
    
    :param pdbString: a string containing a structure in PDB format
    :param considerWaters: take waters into account for energy computations
    
    :return: pandas Dataframe with energy terms of the analyzed structure
             with index corresponding to pair of molecules
    """
    
    if type(pdbString) == type([]):
        pdbString = "\n".join(pdbString)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            df = pd.DataFrame.from_dict(pyFoldXpp.getComplexEnergy(pdbString, considerWaters), orient="index", columns=ENERGY_TERMS)
    else:
            df = pd.DataFrame.from_dict(pyFoldXpp.getComplexEnergy(pdbString, considerWaters), orient="index", columns=ENERGY_TERMS)
    
    return df

def getResiduesEnergy(pdbString, considerWaters=False):
    """
    Compute the interface energy of FoldX for molecules within a structure
    
    :param pdbString: a string containing a structure in PDB format
    :param considerWaters: take waters into account for energy computations
    
    :return: pandas Dataframe with energy terms of the analyzed structure
             with index corresponding to each residue
    """
    
    if type(pdbString) == type([]):
        pdbString = "\n".join(pdbString)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            df = pd.DataFrame.from_dict(pyFoldXpp.getResiduesEnergy(pdbString, considerWaters), orient="index", columns=ENERGY_TERMS)
    else:
        df = pd.DataFrame.from_dict(pyFoldXpp.getResiduesEnergy(pdbString, considerWaters), orient="index", columns=ENERGY_TERMS)
    
    return df

def mutate(pdbString, mutations, numberOfRuns=1):
    """
    Generate a mutated structure
    
    :param pdbString: a string containing a structure in PDB format
    :param mutations: mutations to be performed in foldx format (example GI1A means to mutate GLY in position 1 of molecule I to ALA)
    
    :return: tuple, where element 0 is the DataFrame of ddGs of each models respect to wildtype structure 
            and element 1 is the array of strings with the mutated model
    """
    
    if type(pdbString) == type([]):
        pdbString = "\n".join(pdbString)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            ddGs, mutModels = pyFoldXpp.getMutant(pdbString, mutations, numberOfRuns) 
    else:
        ddGs, mutModels = pyFoldXpp.getMutant(pdbString, mutations, numberOfRuns) 
    
    return (pd.DataFrame(ddGs, columns=ENERGY_TERMS), mutModels)
    
def repair(pdbString,fixResidues=[]):
    """
    Repair the sidechains of a structure
    
    :param pdbString: a string containing a structure in PDB format
    :param fixResidues: list of residues to remain fixed in foldx format (example GI1 means GLY in position 1 of molecule I)
    
    :return: repaired model in string format
    """
    
    if type(pdbString) == type([]):
        pdbString = "\n".join(pdbString)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            ret = pyFoldXpp.getRepairedStructure(pdbString, fixResidues)
    else:
        ret = pyFoldXpp.getRepairedStructure(pdbString, fixResidues)
    
    return ret

def getNetworks(pdbString):
    
    """
    Get networking information of a structure
    
    :param pdbString: a string containing a structure in PDB format

    :return: dictionary where keys are network types and elements are linked residue names for that network
    """
    
    if type(pdbString) == type([]):
        pdbString = "\n".join(pdbString)
    
    if not in_notebook():
        out = OutputGrabber()
        with out: 
            ret = pyFoldXpp.getNetworks(pdbString)
    else:
        ret = pyFoldXpp.getNetworks(pdbString)
    
    return ret

if __name__ == "__main__":
    code = "2ci2"
    #st = Structure(code)
    #print( getTotalEnergy(code, st.toPdb(), False) )
    #print( getResiduesEnergy(st.toPdb(), False) )
    #print( mutate(st.toPdb(), "GI29A;",1)[0] )
    #print( mutate(st.toPdb(), "GI29A;",1)[1][0][0] )
    
    #print( repair(st.toPdb(), ["GI29"]) )
    #for k,v in getNetworks(st.toPdb()).items() : print( k , v )
    '''
    code = "5XJL"
    st = Structure(code)
    print (getInterfaceEnergy(st.toPdb(), False) )
    '''
    
    
    
    
    