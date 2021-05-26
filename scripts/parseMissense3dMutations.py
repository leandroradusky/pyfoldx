'''
Created on Apr 21, 2021
@author: lradusky
'''

import sys
sys.path.append("/home/lradusky/Dropbox/workspacesbg/pyFoldX")

from pyfoldx.structure.structure import Structure
from pyfoldx.handlers.fileHandler import FileHandler
from pyfoldx.structure.misc import OneThree
from pyfoldx.handlers.systemHandler import SystemHandler

from glob import glob
import pandas as pd
import warnings
import numpy as np
from multiprocessing import Pool

warnings.filterwarnings('ignore')

SOURCE_PATH = "/users/lserrano/lradusky/pdb-missense/"
DEST_PATH = "/users/lserrano/lradusky/Missense-processed/"

def getDDG(mutant_prefix):
    ddG_df = pd.read_csv(mutant_prefix+"_mutDiff.csv", header=0, index_col=0)
    ddG_df.columns = ["#DDG_"+ str(col) for col in ddG_df.columns]
    return ddG_df

def getResidueDetails(mutant_prefix, mutant):
    
    mt_resDet_df = pd.read_csv(mutant_prefix+"_residuesEnergy.csv", header=0, index_col=None)
    
    ori_res, chain, pos, mut_res = mutant[0], mutant[1], mutant[2:-1], mutant[-1]
    ret_df = mt_resDet_df[mt_resDet_df["Pos"] == int(pos)][mt_resDet_df["Mol"] == chain][["omega","phi","psi","sec_struct"]]
    ret_df.columns = ["#RES_"+ str(col) for col in ret_df.columns]
    
    return ret_df

def getInterfaceDetails(wt_prefix, mutant_prefix, mutant, pdb):
    
    ori_res, chain, pos, mut_res = mutant[0], mutant[1], mutant[2:-1], mutant[-1]
    
    if not FileHandler.fileExists(wt_prefix+"_interfaceEnergy.csv"): return 0, 0, 0
    
    mt_res_in = int( FileHandler.getLines(mutant_prefix+"_interfaceResidues.obj")[0].find("""'%s%s%s'""" % (mut_res, chain, pos) ) != -1 )
    wt_res_in = int( FileHandler.getLines(wt_prefix+"_interfaceResidues.obj")[0].find("""'%s%s%s'""" % (ori_res, chain, pos) ) != -1 )
    
    if mt_res_in == 0 and wt_res_in == 0: return 0, 0, 0
    
    st_wt = Structure("WT", SOURCE_PATH+"T_/WT_"+pdb[0]+"/WT_"+pdb+"_"+mutant+"_0_WT.pdb")
    st_wt.getInterfaceEnergy(verbose=False)

    mt_if_energy_df =  pd.read_csv(mutant_prefix+"_interfaceEnergy.csv", header=0, index_col=[0,1])
    
    st_if_energy_df = st_wt.interfaceEnergy[["Interaction Energy"]] 
    mt_if_energy_df = mt_if_energy_df[["Interaction Energy"]]
    
    st_if_energy_df = st_if_energy_df.iloc[st_if_energy_df.index.get_level_values('Group1') == chain].append(\
                      st_if_energy_df.iloc[st_if_energy_df.index.get_level_values('Group2') == chain])
    mt_if_energy_df = mt_if_energy_df.iloc[mt_if_energy_df.index.get_level_values('Group1') == chain].append(\
                      mt_if_energy_df.iloc[mt_if_energy_df.index.get_level_values('Group2') == chain])
    
    int_ddG = float((mt_if_energy_df.astype(float) - st_if_energy_df.astype(float)).sum()) 
    
    return int_ddG , mt_res_in, wt_res_in 

def getNetworkDetails(wt_prefix, mutant_prefix, mutant, pdb):
    
    ori_res, chain, pos, mut_res = mutant[0], mutant[1], mutant[2:-1], mutant[-1]
    
    st_wt = Structure("WT", SOURCE_PATH+"T_/WT_"+pdb[0]+"/WT_"+pdb+"_"+mutant+"_0_WT.pdb")
    st_wt.getNetworks(verbose=False)
    
    mt_vdw_df = pd.read_csv(mutant_prefix+"_networks_VdWClashes.csv", header=0, index_col=0)
    wt_vdw_df = st_wt.networks["VdWClashes"]
    mt_vdw_df.index = [x[3:] for x in mt_vdw_df.index ]
    wt_vdw_df.index = [x[3:] for x in wt_vdw_df.index ]
    vdwClashes = mt_vdw_df.loc[chain+pos].sum() 
    vdwClashes_diff = mt_vdw_df.loc[chain+pos].sum() - wt_vdw_df.loc[chain+pos].sum()
    
    mt_hbo_df = pd.read_csv(mutant_prefix+"_networks_Hbonds.csv", header=0, index_col=0)
    wt_hbo_df = st_wt.networks["Hbonds"]
    mt_hbo_df.index = [x[3:] for x in mt_hbo_df.index ]
    wt_hbo_df.index = [x[3:] for x in wt_hbo_df.index ]
    Hbonds = mt_hbo_df.loc[chain+pos].sum() 
    Hbonds_diff = mt_hbo_df.loc[chain+pos].sum() - wt_hbo_df.loc[chain+pos].sum()
    
    return vdwClashes, vdwClashes_diff, Hbonds, Hbonds_diff
    
def computeRelativeBFactor(wt_prefix, mutant, pdb):
    
    try:
        ori_res, chain, pos, mut_res = mutant[0], mutant[1], mutant[2:-1], mutant[-1]
        
        st = Structure(pdb, path = wt_prefix + ".pdb")
        
        min_bf, max_bf = st.getMinMaxBFactor()
        avg_res_bf = st.data[chain][int(pos)].getAvgBFactor()
        
        return (float(avg_res_bf) - float(min_bf)) / (float(max_bf) - float(min_bf))
    except:
        return np.nan

def computeWaterMediated():
    pass

def processPdb(pdb):
    
    try:
        working_path = SOURCE_PATH+pdb[1:3]+"/"+pdb+"/"
        wt_st_prefix = working_path+pdb
        
        ### Create the  dataframe of this pdb ###
        pdb_df = pd.DataFrame(columns=['#DDG_total', '#DDG_backHbond', '#DDG_sideHbond', '#DDG_energy_VdW', 
                                       '#DDG_electro', '#DDG_energy_SolvP', '#DDG_energy_SolvH', '#DDG_energy_vdwclash', 
                                       '#DDG_entrop_sc', '#DDG_entrop_mc', '#DDG_sloop_entropy', '#DDG_mloop_entropy', 
                                       '#DDG_cis_bond', '#DDG_energy_torsion', '#DDG_backbone_vdwclash', '#DDG_energy_dipole', 
                                       '#DDG_water', '#DDG_disulfide', '#DDG_energy_kon', '#DDG_partcov', '#DDG_energyIonisation', 
                                       '#DDG_entr_complex', '#RES_omega', '#RES_phi', '#RES_psi', '#RES_sec_struct', '#INTERFACE_ddG', 
                                       '#INTERFACE_wt_res', '#INTERFACE_mt_res', '#NETWORKS_vdw', '#NETWORKS_vdw_diff', '#NETWORKS_hbond', 
                                       '#NETWORKS_hbond_diff', '#BFACTOR_relative'])
        
        for m in glob( working_path+pdb+"_*.pdb" ):
            
            mutant = m.split("_")[-2]
            
            print( "***********" )
            print( pdb, mutant )
            
            
            try:
                mutant_st_prefix = working_path+pdb+"_"+mutant+"_0" 
                
                ### Copy mutant files to new folder ###
                
                working_destination = DEST_PATH+pdb[1:3]+"/"+pdb+"/"
                FileHandler.ensureDir(working_destination)
                SystemHandler.executeCommand("cp %s %s" % (mutant_st_prefix+".pdb", working_destination+pdb+"_"+mutant+".pdb"))
                
                ### Parse ddG ###
                
                ddG_df = getDDG(mutant_st_prefix)
                #print ( ddG_df )
                
                ### Parse residue details ###
                
                res_df = getResidueDetails(mutant_st_prefix, mutant)
                #print ( res_df )
                
                ### Parse interface details ###
                
                int_ddG , mt_res_is_interface, wt_res_is_interface  = getInterfaceDetails(wt_st_prefix, mutant_st_prefix, mutant, pdb)
                #print ( int_ddG , mt_res_is_interface, wt_res_is_interface )
                
                ### Parse network details ###
                
                vdwClashes, vdwClashes_diff, Hbonds, Hbonds_diff = getNetworkDetails(wt_st_prefix, mutant_st_prefix, mutant, pdb)
                #print( vdwClashes, vdwClashes_diff, Hbonds, Hbonds_diff, float(ddG_df["#DDG_total"]) )
                
                ### Compute relative beta factor ###
                
                rel_bfactor = computeRelativeBFactor(wt_st_prefix, mutant, pdb)
                #print( rel_bfactor )
                
                ### Compute if water mediated ###
                
                ### Generate JSon for PDBe ###
                
                ### Join everything to make a row for this mutation ###
                resDict= ddG_df.to_dict(orient="list") 
                resDict.update( res_df.to_dict(orient="list") )
                resDict.update( {"#INTERFACE_ddG": int_ddG, "#INTERFACE_wt_res": wt_res_is_interface, "#INTERFACE_mt_res": mt_res_is_interface} ) 
                resDict.update( {"#NETWORKS_vdw": vdwClashes, "#NETWORKS_vdw_diff": vdwClashes_diff} ) 
                resDict.update( {"#NETWORKS_hbond": Hbonds, "#NETWORKS_hbond_diff": Hbonds_diff} ) 
                resDict.update( {"#BFACTOR_relative": rel_bfactor} ) 
                
                ### Append it into the dataframe of this pdb ###
                pdb_df.loc[mutant] = pd.DataFrame.from_dict( resDict ).loc[0] 
            
            except:
                print( "Failed" )
            
            print( "***********" )
        
        ### Save the  dataframe of this pdb ###
        pdb_df.to_csv(DEST_PATH+"CSV/"+pdb+".csv")
        FileHandler.appendLine(DEST_PATH+"processed", pdb)
    except:
        FileHandler.appendLine(DEST_PATH+"failed", pdb)


def generatePdbCsv():
    
    processed = set(FileHandler.getLines(SOURCE_PATH+"processed"))
    processed = processed - set(FileHandler.getLines(DEST_PATH+"processed")) 
    processed = sorted(processed)
    
    print ("Processing %s pdbs" % len(processed))
    
    # Parallel version
    pool = Pool(processes=6)
    pool.map(processPdb, processed)
    
    # Serial version
    #for pdb in processed:
    #    processPdb(pdb)
    
    print(len(processed))
    
def generateJointDataset():
    
    miss_3d_path = "/home/lradusky/Dropbox/raduspostdoc/pyFoldXpaper/DiseaseMutations/missense3d-benchmarking_PDBeKB.csv"
    miss_3d_dest = "/home/lradusky/Dropbox/raduspostdoc/pyFoldXpaper/DiseaseMutations/missense3d-benchmarking_PDBeKB_foldx.csv"
    
    ms3d_df = pd.read_csv(miss_3d_path, header=0)
    ms3d_new = pd.read_csv(miss_3d_path, header=0)
    
    d = False

    for i, row in ms3d_df.iterrows():
        pdb = row["#PDB"]
        if not FileHandler.fileExists(DEST_PATH+"CSV/"+pdb+".csv"): continue
        
        pdb_df = pd.read_csv(DEST_PATH+"CSV/"+pdb+".csv", header=0, index_col=0)
        
        mutation = row["#RESWT"]+row["#CHAIN"]+str(row["#PDBPOS"])+row["#RESMUT"]
        
        if not d:
            for col in pdb_df.columns:
                ms3d_new.loc[:, col] = np.nan
            d=True
        
        try: ms3d_new.loc[i, pdb_df.columns] = pdb_df.loc[mutation]
        except: pass
        
        print(i)
        
    ms3d_new.to_csv(miss_3d_dest)

if __name__ == "__main__":
    
    print( "started" )
    
    #generatePdbCsv()
    generateJointDataset()
    
    print( "done" )
    