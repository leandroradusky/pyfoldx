'''
Created on Mar 18, 2021
@author: lradusky
'''

import pandas as pd

from pyfoldx.handlers.fileHandler import FileHandler
from pyfoldx.handlers.systemHandler import SystemHandler
import pyfoldx.structure.structure as structure
from pyfoldx.structure.misc import in_notebook, OutputGrabber
from _collections import defaultdict
from datetime import datetime
import os

ENERGY_TERMS = ['total','backHbond','sideHbond','energy_VdW','electro','energy_SolvP','energy_SolvH','energy_vdwclash',
                'entrop_sc','entrop_mc','sloop_entropy','mloop_entropy','cis_bond','energy_torsion','backbone_vdwclash',
                'energy_dipole','water','disulfide','energy_kon','partcov','energyIonisation','entr_complex']

FOLDX_LOCATION = os.getenv('FOLDX_BINARY')

def getJobFolder():
    '''
    Generates a folder with a name based on a timestamp to allow unique folders for paralellization
    '''
    a = datetime.utcnow()
    jobid=str(a.year)+str(a.month)+str(a.day)+str(a.hour)+str(a.minute)+str(a.second)+str(a.microsecond)
    return "./.foldx_%s/" % jobid

def checkFoldx():
    if not FileHandler.fileExists(FOLDX_LOCATION):
        print("FoldX executable not found, please set the variable FOLDX_LOCATION")
        exit(0)
    return True

def runFoldX(command,working_folder, pdb,other_parameters={}, working_with_rna=False, silent=True, consider_waters=False):
    
    if FileHandler.exists("./molecules/"):
        SystemHandler.executeCommand("cp -r molecules %s" % working_folder)
        
    if not in_notebook():
        out = OutputGrabber()
    else:
        out = open(working_folder+"./tmp.pdb","r")
    
    with out: 
        checkFoldx()
        
        
        c="cd "+working_folder+"; "
        c+=FOLDX_LOCATION+" --command="+command+" --pdb="+pdb
        
        if consider_waters:
            other_parameters["water"] = "-CRYSTAL"
            other_parameters["pdbWaters"]="true" 
        
        for paramName, paramValue in other_parameters.items():
            c+= " --%s=%s " % (paramName, paramValue)
        if working_with_rna:
            c+=" --complexWithRNA=true "
        if silent:
            c+=" > /dev/null 2> /dev/null"
        if not silent:
            print( "Running command: " )
            print( c )
        SystemHandler.executeCommand(c)

def getTotalEnergy(st, consider_waters=False, other_parameters={}):
    """
    Compute the stability energy of FoldX for a structure
    
    :param st: structure object to analyze
    :param consider_waters: take waters into account for energy computations
    :param other_parameters: special foldx parameters
    
    :return: pandas Dataframe with energy terms of the analyzed structure
    
    :todo: unrecognized molecules
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    # Run FoldX
    runFoldX(command          = "Stability", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = consider_waters )
    
    # Read files
    st_lines = FileHandler.getLines(TMP_FOLDER+"tmp_0_ST.fxout")
    Unrecognized_molecules = FileHandler.getLines(TMP_FOLDER+"Unrecognized_molecules.txt")
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    # Parse files
    dfData = st_lines[0].split("\t")[1:-1]
    st.totalEnergy= pd.DataFrame(data = [dfData], columns= ENERGY_TERMS, index=[st.code])
    
    # Return Dataframe
    return st.totalEnergy

def getComplexEnergy(st, consider_waters=False, other_parameters={}):
    """
    Compute the molecule pairwise interaction energy of FoldX for a structure
    
    :param st: structure object to analyze
    :param consider_waters: take waters into account for energy computations
    :param other_parameters: special foldx parameters
    
    :return: pandas Dataframe with the pairwise interaction energy terms within the analyzed structure.
    
    :todo: unrecognized molecules and the rest of the dfs
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    runFoldX(command          = "AnalyseComplex", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = consider_waters )
    
    # Read files
    Indiv_energies = FileHandler.getLines(TMP_FOLDER+"Indiv_energies_tmp_AC.fxout")
    Interface_Residues = FileHandler.getLines(TMP_FOLDER+"Interface_Residues_tmp_AC.fxout")
    Interaction = FileHandler.getLines(TMP_FOLDER+"Interaction_tmp_AC.fxout")
    #Summary = fileHandler.getLines(TMP_FOLDER+"Summary_tmp_AC.fxout")
    Unrecognized_molecules = FileHandler.getLines(TMP_FOLDER+"Unrecognized_molecules.txt")
    
    #Parse files
    interactionDfData = [ tuple(x.split("\t")[1:]) for x in Interaction[9:] ]
    interactionDfColumns = Interaction[8].split("\t")[1:]
    st.interfaceEnergy = pd.DataFrame(data = interactionDfData, columns= interactionDfColumns)
    st.interfaceEnergy.set_index(["Group1","Group2"], inplace=True)
    
    IndivDfData = [ tuple(x.split("\t")[1:]) for x in Indiv_energies[9:] ]
    IndivDfColumns = Indiv_energies[8].split("\t")[1:]
    st.molEnergy = pd.DataFrame(data = IndivDfData, columns= IndivDfColumns)
    st.molEnergy.set_index(["Group"], inplace=True)
    
    st.interfaceResidues = {}
    start=False
    for i in range(len(Interface_Residues)):
        if Interface_Residues[i].find("interface residues between") != -1:
            start=True
            mols = (Interface_Residues[i].split()[-3],Interface_Residues[i].split()[-1])
        elif start==True:
            st.interfaceResidues[mols]= Interface_Residues[i].split()
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    return st.interfaceEnergy

def getResiduesEnergy(st, consider_waters, other_parameters={}):
    """
    Compute the per residue detailed energy of FoldX for a structure
    
    :param st: structure object to analyze
    :param consider_waters: take waters into account for energy computations
    :param other_parameters: special foldx parameters
    
    :return: pandas Dataframe with the pairwise interaction energy terms within the analyzed structure.
    
    :todo: unrecognized molecules and the rest of the dfs
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    # Run FoldX
    runFoldX(command          = "SequenceDetail", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = consider_waters )
    
    # Read files
    sd_lines = FileHandler.getLines(TMP_FOLDER+"SD_tmp.fxout")
    Unrecognized_molecules = FileHandler.getLines(TMP_FOLDER+"Unrecognized_molecules.txt")
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    # Parse files
    dfData = { line.split("\t")[-1] : line.split("\t")[1:-7] for line in sd_lines }
    st.residuesEnergy = pd.DataFrame.from_dict(dfData, orient="index", columns=["Code", "Mol", "Pos","omega","phi","psi","sec_struct"]+ENERGY_TERMS)
    st.residuesEnergy.set_index(["Code", "Mol", "Pos"], inplace=True)
    st.residuesEnergy.dropna(axis = 0, how = 'all', inplace = True)
    
    # Return Dataframe
    return st.residuesEnergy
    
def getMutants(st, mutations, number_of_runs=1, other_parameters={}):
    """
    Generate a mutated structure(s)
    
    :param st: structure object to mutate
    :param mutations: mutations to be performed in foldx format (example GI1A means to mutate GLY in position 1 of molecule I to ALA)
    :param number_of_runs: number of mutations to generate
    :param other_parameters: special foldx parameters
    
    :return: tuple, where element 0 is the DataFrame of ddGs of each models respect to wildtype structure 
            element 1 is the list of mutated structures and element 2 is the list of wt structures.
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    # Write individual_list
    FileHandler.writeLine(TMP_FOLDER+"individual_list.txt", mutations)
    
    # Set parameters of BuildModel
    other_parameters["numberOfRuns"] = str(number_of_runs)
    other_parameters["mutant-file"] = "individual_list.txt"
    
    # Run FoldX
    runFoldX(command          = "BuildModel", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = False )
    
    # Read Files
    mutPdbs= []
    wtPdbs= []
    for i in range(number_of_runs):
        mutPdbs.append( structure.Structure(st.code+"_%s_%s" % (mutations.replace(",","_").replace(";",""),i), 
                                            TMP_FOLDER+"tmp_1%s.pdb" % ("_"+str(i) if number_of_runs > 1 else "") ) )
        wtPdbs.append( structure.Structure("WT_"+st.code+"_%s_%s" % (mutations.replace(",","_").replace(";",""), i),
                                            TMP_FOLDER+"WT_tmp_1%s.pdb" % ("_"+str(i) if number_of_runs > 1 else "") ) )
    
    raw_lines = FileHandler.getLines(TMP_FOLDER+"Raw_tmp.fxout")
    for i in range(number_of_runs):
        mut_line_number = number_of_runs*(-2) + i
        dfData = raw_lines[mut_line_number].split("\t")[1:]
        mutPdbs[i].totalEnergy = pd.DataFrame(data = [dfData], columns= ENERGY_TERMS, index=[mutPdbs[i].code])
        
        wt_line_number = -number_of_runs + i
        dfData = raw_lines[wt_line_number].split("\t")[1:]
        wtPdbs[i].totalEnergy = pd.DataFrame(data = [dfData], columns= ENERGY_TERMS, index=[wtPdbs[i].code])
    
    dif_lines = FileHandler.getLines(TMP_FOLDER+"Dif_tmp.fxout")
    difData = {}
    for i in range(number_of_runs):
        mut_line_number = -number_of_runs + i
        difData[st.code+"_%s_%s" % (mutations.replace(",","_").replace(";",""), i)] = dif_lines[mut_line_number].split("\t")[1:]
    
    retDf = pd.DataFrame.from_dict(data = difData, orient="index", columns= ENERGY_TERMS)
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    return retDf, mutPdbs, wtPdbs

def getRepairedStructure(st, fix_residues=[], other_parameters={}):
    """
    Repair the sidechains of a structure
    
    :param pdb_string: a string containing a structure in PDB format
    :param fix_residues: list of residues to remain fixed in foldx format (example GI1 means GLY in position 1 of molecule I)
    
    :return: repaired model in string format
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    # Set parameters of BuildModel
    if fix_residues != []:
        FileHandler.writeLine(TMP_FOLDER+"fix-residues.txt", ",".join(fix_residues))
        other_parameters["fix-residues-file"] = "fix-residues.txt"
    
    # Run FoldX
    runFoldX(command          = "RepairPDB", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = False )
    
    # Read Files
    repSt = structure.Structure(st.code+"_Repaired", TMP_FOLDER+"tmp_Repair.pdb" )
    Unrecognized_molecules = FileHandler.getLines(TMP_FOLDER+"Unrecognized_molecules.txt")
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    return repSt

def getNetworks(st, other_parameters={}):
    """
    Get the networks of a structure
    
    :param st: structure object to get networks from
    
    :summary: set the network attributes of a structure, no return.
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    # Run FoldX
    runFoldX(command          = "PrintNetworks", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = False )
    
    # Read Files
    st.networks = {}
    n_residues = len(FileHandler.getLines(TMP_FOLDER+"Matrix_Distances_tmp_PN.fxout")[4].split("\t"))-1
    filfunc = lambda x:x<4 or x>4+n_residues
    st.networks["Distances"]  = pd.read_csv(TMP_FOLDER+"Matrix_Distances_tmp_PN.fxout" , header=0, index_col=0, skiprows=filfunc, sep="\t")
    st.networks["Disulfide"]  = pd.read_csv(TMP_FOLDER+"Matrix_Disulfide_tmp_PN.fxout" , header=0, index_col=0, skiprows=filfunc, sep="\t")
    st.networks["Electro"]    = pd.read_csv(TMP_FOLDER+"Matrix_Electro_tmp_PN.fxout"   , header=0, index_col=0, skiprows=filfunc, sep="\t")
    st.networks["Hbonds"]     = pd.read_csv(TMP_FOLDER+"Matrix_Hbonds_tmp_PN.fxout"    , header=0, index_col=0, skiprows=filfunc, sep="\t")
    st.networks["Partcov"]    = pd.read_csv(TMP_FOLDER+"Matrix_Partcov_tmp_PN.fxout"   , header=0, index_col=0, skiprows=filfunc, sep="\t")
    st.networks["VdWClashes"] = pd.read_csv(TMP_FOLDER+"Matrix_VdWClashes_tmp_PN.fxout", header=0, index_col=0, skiprows=filfunc, sep="\t")
    st.networks["Volumetric"] = pd.read_csv(TMP_FOLDER+"Matrix_Volumetric_tmp_PN.fxout", header=0, index_col=0, skiprows=filfunc, sep="\t")
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
def getAlaScan(st, consider_waters, other_parameters={}):
    """
    Compute the energy variation upon mutation to alanine of FoldX for a structure
    
    :param st: structure object to analyze
    :param consider_waters: take waters into account for energy computations
    :param other_parameters: special foldx parameters
    
    :return: pandas Dataframe with the pairwise interaction energy terms within the analyzed structure.
    
    :todo: unrecognized molecules and the rest of the dfs
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    # Run FoldX
    runFoldX(command          = "AlaScan", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = consider_waters )
    
    # Read Files
    as_lines = FileHandler.getLines(TMP_FOLDER+"tmp_AS.fxout")
    as_data = [(x.split()[0]+"_"+x.split()[1],x.split()[-1] ) for x in as_lines]
    ret = pd.DataFrame(data=as_data, columns=["Residue","ddG_ala"])
    ret.set_index(["Residue"],inplace=True)
    
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    return ret

def getPSSM(st, positions_string, other_parameters={}):
    """
    Compute the position specific scoring matrix (PSSM) of FoldX for a structure
    
    :param st: structure object to analyze
    :param other_parameters: special foldx parameters
    
    :return: pandas Dataframe with the pairwise interaction energy terms within the analyzed structure.
    
    :todo: unrecognized molecules and the rest of the dfs
    """
    # Get folder name for this job
    TMP_FOLDER = getJobFolder()
    
    # Write pdb
    st.toPdbFile(TMP_FOLDER+"tmp.pdb")
    
    other_parameters["positions"] = positions_string
    other_parameters["out-pdb"] = "false"
    
    # Run FoldX
    runFoldX(command          = "Pssm", 
             working_folder   = TMP_FOLDER, 
             pdb              = "tmp.pdb",
             other_parameters = other_parameters,
             consider_waters  = False )
    
    # Read Files
    mutsLines = FileHandler.getLines(TMP_FOLDER+"individual_list_0_PSSM.txt")
    energies = pd.read_csv(TMP_FOLDER+"Dif__0_tmp.fxout" , header=0, index_col=0, skiprows=8, sep="\t")
    ret = defaultdict(lambda:defaultdict(lambda:"Nan"))
    
    mutN = 1
    for line in mutsLines:
        resCode = line[0]
        resChain = line[1]
        resNum = line[2:-2]
        resMut = line[-2]
        
        energy = energies.loc["tmp_%s.pdb" % mutN]["total energy"]
        
        ret[resCode+resNum][resMut] = energy
        
        mutN+=1
    
    ret = pd.DataFrame.from_dict(data=ret)
    # Delete files
    SystemHandler.executeCommand("rm -r %s" % TMP_FOLDER)
    
    return ret
    
    
if __name__ == "__main__":
    pass



