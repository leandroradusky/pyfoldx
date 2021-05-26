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

SOURCE_PATH = "/users/lserrano/lradusky/frustration/"

def getFrustration(pdb, mutant, mutant_path, resolution):
    
    ori_res, chain, pos, mut_res = mutant[0], mutant[1], mutant[2:-1], mutant[-1]
    wt_path = mutant_path.replace("_"+mutant+"_0","")
    
    wt_confi_path = wt_path+"/FrustrationData/"+pdb+".pdb_configurational_5adens"
    wt_mutat_path = wt_path+"/FrustrationData/"+pdb+".pdb_mutational_5adens"
    wt_singl_path = wt_path+"/FrustrationData/"+pdb+".pdb_singleresidue"
    
    mut_confi_path = mutant_path+"/FrustrationData/"+pdb+"_"+mutant+"_0.pdb_configurational_5adens"
    mut_mutat_path = mutant_path+"/FrustrationData/"+pdb+"_"+mutant+"_0.pdb_mutational_5adens"
    mut_singl_path = mutant_path+"/FrustrationData/"+pdb+"_"+mutant+"_0.pdb_singleresidue"
    
    if FileHandler.fileExists(wt_confi_path) and \
       FileHandler.fileExists(wt_mutat_path) and \
       FileHandler.fileExists(wt_singl_path) and \
       FileHandler.fileExists(mut_confi_path) and \
       FileHandler.fileExists(mut_mutat_path) and \
       FileHandler.fileExists(mut_singl_path):
        
        # Single Residue
        command= 'grep "^%s %s" %s' % (pos, chain, wt_singl_path)
        sr_wt_Res, sr_wt_ChainRes, sr_wt_DensityRes, sr_wt_AA, sr_wt_NativeEnergy, sr_wt_DecoyEnergy, sr_wt_SDEnergy, sr_wt_FrstIndex = \
            SystemHandler.getCommandResult(command)[0].split()
        command= 'grep "^%s %s" %s' % (pos, chain, mut_singl_path)
        sr_mt_Res, sr_mt_ChainRes, sr_mt_DensityRes, sr_mt_AA, sr_mt_NativeEnergy, sr_mt_DecoyEnergy, sr_mt_SDEnergy, sr_mt_FrstIndex = \
            SystemHandler.getCommandResult(command)[0].split()
        
        # Configurational
        command= 'grep "^%s %s" %s' % (pos, chain, wt_confi_path)
        confi_wt_Res, confi_wt_ChainRes, confi_wt_Total, confi_wt_nHighlyFrst, confi_wt_nNeutrallyFrst, \
        confi_wt_nMinimallyFrst, confi_wt_relHighlyFrustrated, confi_wt_relNeutralFrustrated, \
        confi_wt_relMinimallyFrustrated = SystemHandler.getCommandResult(command)[0].split()
        
        command= 'grep "^%s %s" %s' % (pos, chain, mut_confi_path)
        confi_mt_Res, confi_mt_ChainRes, confi_mt_Total, confi_mt_nHighlyFrst, confi_mt_nNeutrallyFrst, \
        confi_mt_nMinimallyFrst, confi_mt_relHighlyFrustrated, confi_mt_relNeutralFrustrated, \
        confi_mt_relMinimallyFrustrated = SystemHandler.getCommandResult(command)[0].split()
        
        # Mutational
        command= 'grep "^%s %s" %s' % (pos, chain, wt_mutat_path)
        mutat_wt_Res, mutat_wt_ChainRes, mutat_wt_Total, mutat_wt_nHighlyFrst, mutat_wt_nNeutrallyFrst, \
        mutat_wt_nMinimallyFrst, mutat_wt_relHighlyFrustrated, mutat_wt_relNeutralFrustrated, \
        mutat_wt_relMinimallyFrustrated = SystemHandler.getCommandResult(command)[0].split()
        
        command= 'grep "^%s %s" %s' % (pos, chain, mut_mutat_path)
        mutat_mt_Res, mutat_mt_ChainRes, mutat_mt_Total, mutat_mt_nHighlyFrst, mutat_mt_nNeutrallyFrst, \
        mutat_mt_nMinimallyFrst, mutat_mt_relHighlyFrustrated, mutat_mt_relNeutralFrustrated, \
        mutat_mt_relMinimallyFrustrated = SystemHandler.getCommandResult(command)[0].split()
        
        line = "\t".join([ str(x) for x in [pdb, mutant, resolution, \
        sr_wt_Res, sr_wt_ChainRes, sr_wt_DensityRes, sr_wt_AA, sr_wt_NativeEnergy, sr_wt_DecoyEnergy, sr_wt_SDEnergy, sr_wt_FrstIndex,\
        sr_mt_Res, sr_mt_ChainRes, sr_mt_DensityRes, sr_mt_AA, sr_mt_NativeEnergy, sr_mt_DecoyEnergy, sr_mt_SDEnergy, sr_mt_FrstIndex,\
        confi_wt_Res, confi_wt_ChainRes, confi_wt_Total, confi_wt_nHighlyFrst, confi_wt_nNeutrallyFrst, \
        confi_wt_nMinimallyFrst, confi_wt_relHighlyFrustrated, confi_wt_relNeutralFrustrated, confi_wt_relMinimallyFrustrated, \
        confi_mt_Res, confi_mt_ChainRes, confi_mt_Total, confi_mt_nHighlyFrst, confi_mt_nNeutrallyFrst, \
        confi_mt_nMinimallyFrst, confi_mt_relHighlyFrustrated, confi_mt_relNeutralFrustrated, confi_mt_relMinimallyFrustrated, \
        mutat_wt_Res, mutat_wt_ChainRes, mutat_wt_Total, mutat_wt_nHighlyFrst, mutat_wt_nNeutrallyFrst, \
        mutat_wt_nMinimallyFrst, mutat_wt_relHighlyFrustrated, mutat_wt_relNeutralFrustrated, mutat_wt_relMinimallyFrustrated, \
        mutat_mt_Res, mutat_mt_ChainRes, mutat_mt_Total, mutat_mt_nHighlyFrst, mutat_mt_nNeutrallyFrst, \
        mutat_mt_nMinimallyFrst, mutat_mt_relHighlyFrustrated, mutat_mt_relNeutralFrustrated, mutat_mt_relMinimallyFrustrated ] ] )
        
        FileHandler.appendLine(SOURCE_PATH+"frustration.tsv", line)
        
    return

def processPdb(pdb):
    
    try:
        
        ### Get Structure Resolution ###
        
        resolution = Structure(pdb).resolution
        
        working_path = SOURCE_PATH+pdb[1:3]+"/"+pdb+"/"
        wt_st_prefix = working_path+pdb
        
        for m in glob( working_path+pdb+"_*.done" ):
            
            mutant = m.split("_")[-2]
            
            #print( "***********" )
            print( pdb, mutant, resolution )
            
            ### Get Frustration values ###
            getFrustration(pdb, mutant, m, resolution)
            
            
            #print( "***********" )
        
        FileHandler.appendLine(SOURCE_PATH+"parsed", pdb)
    except:
        FileHandler.appendLine(SOURCE_PATH+"parsefailed", pdb)


def generateParsedFile():
    
    processed = set(FileHandler.getLines(SOURCE_PATH+"processed"))
    processed = processed - set(FileHandler.getLines(SOURCE_PATH+"parsed")) 
    processed = sorted(processed)
    
    print ("Processing %s pdbs" % len(processed))
    
    # Parallel version
    #pool = Pool(processes=6)
    #pool.map(processPdb, processed)
    
    # Serial version
    j=0
    for pdb in processed:
        processPdb(pdb)
        j+=1
        print ("processed %s of %s" % (j, len(processed)))
    
    print(len(processed))

if __name__ == "__main__":
    
    print( "started" )
    
    generateParsedFile()
    
    print( "done" )
    