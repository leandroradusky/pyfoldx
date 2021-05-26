'''
Created on Apr 9, 2021
@author: lradusky
'''

import sys
sys.path.append("/home/lradusky/Dropbox/workspacesbg/pyFoldX")

from pyfoldx.structure.structure import Structure
from pyfoldx.handlers.fileHandler import FileHandler
import pandas as pd
from multiprocessing import Pool

MISSENSE_3D_PATH = "/home/lradusky/Dropbox/raduspostdoc/pyFoldXpaper/DiseaseMutations/missense3d-benchmarking_PDBeKB.csv"
DESTINATION_PATH = "/users/lserrano/lradusky/pdb-missense/"
#DESTINATION_PATH = "/data/pdb-missense/"
ms3d_df = pd.read_csv(MISSENSE_3D_PATH, header=0)

def saveInfo(st, difDf = None):
    '''
    Given a structure object computes all its foldx info and save it to file
    '''
    
    working_path = DESTINATION_PATH+st.code[1:3]+"/"+st.code[0:4]+"/"
    FileHandler.ensureDir(working_path)
    
    st.toPdbFile(working_path+st.code+".pdb")
    
    st.getTotalEnergy()
    st.getResiduesEnergy()
    st.getNetworks()
    
    st.totalEnergy.to_csv    ( working_path+st.code+"_totalEnergy.csv" )
    st.residuesEnergy.to_csv ( working_path+st.code+"_residuesEnergy.csv" )
    for k in st.networks.keys():
        st.networks[k].to_csv( working_path+st.code+"_networks_%s.csv" % k )
    
    # If not interfaces this could fail we put a try
    try:
        st.getInterfaceEnergy()
        st.interfaceEnergy.to_csv( working_path+st.code+"_interfaceEnergy.csv" )
        st.molEnergy.to_csv      ( working_path+st.code+"_molEnergy.csv" )
        FileHandler.writeLine( working_path+st.code+"_interfaceResidues.obj", str(st.interfaceResidues) ) 
    except:
        print("No interfaces")
    
    if difDf is not None:
        difDf.to_csv( working_path+st.code+"_mutDiff.csv" )
    
def processPdb(pdb):
    if pdb in FileHandler.getLines(DESTINATION_PATH+"processed") :
        # or pdb in FileHandler.getLines(DESTINATION_PATH+"failed"): 
        print ( "Pdb %s processed or failed" % pdb ) 
        return
        
    print ( "Processing pdb %s" % pdb ) 
    
    try:
        # We filter each PDB code to paralellize
        pdb_df = ms3d_df.loc[ ms3d_df["#PDB"] == pdb ]
        # We also instantiate the structure object, to compute and mutate
        st = Structure(pdb)
        
        saveInfo(st)
        
        for dummy, row in pdb_df.iterrows():
            # we form the mutation in foldx format.
            mutation = row["#RESWT"]+row["#CHAIN"]+str(row["#PDBPOS"])+row["#RESMUT"]+";"
            # we mutate
            energyMut, ensMut, enstWt = st.mutate(mutation, verbose = False)
            
            stMut = ensMut.getFrame()
            stWT = enstWt.getFrame()
            
            saveInfo(stMut,energyMut)
            stWT.toPdbFile(DESTINATION_PATH+stWT.code[1:3]+"/"+stWT.code[0:4]+"/"+stWT.code+"_WT.pdb")
        
        FileHandler.appendLine(DESTINATION_PATH+"processed", pdb)
    except:
        FileHandler.appendLine(DESTINATION_PATH+"failed", pdb)
    
def computeMutationsStructure():
    
    # Threaded version
    pool = Pool(processes=7)
    pool.map(processPdb,sorted(set(ms3d_df["#PDB"])))
    # Serial version
    #for pdb in sorted(set(ms3d_df["#PDB"])):
    #    processPdb(pdb)
    
if __name__ == "__main__":
    
    print( "started" )
    
    # compute mutation and save files using the specified structure
    computeMutationsStructure()
    #processPdb("6o1g")
    
    # compute mutation and save files using the specified ensemble
    #computeMutationsEnsemble()
    
    print( "done" )
    
