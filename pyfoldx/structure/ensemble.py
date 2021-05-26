'''
Created on Nov 2, 2020

@author: lradusky
'''

from glob import glob
import pandas as pd
from tqdm import tqdm,trange
import sys
from pyfoldx.handlers.fileHandler import FileHandler
import pyfoldx.structure.structure as structure
import pyfoldx.foldx.foldxHandler as foldxHandler 
from pyfoldx.handlers.uniprot import Orf
from pyfoldx.structure.misc import align

##########################
#### ENSEMBLE
##########################

class Ensemble(object):
    '''
    Class representing ensembles of structures (multiple models)
    
    :ivar frames: list of structure objects.
    :ivar codes: list of codes of the structure objects.
    :ivar chains: list of chains of the structure objects.
    :ivar workingPath: temporary path to download and align structures.
    :ivar code: name of the ensemble.
    '''
    
    def __init__(self, code, working_path, ensemble_file=""):
        '''
        Constructor of an ensemble. 
        
        :param code: Code of the ensemble, just for naming purposes.
        :param working_path: Temporary path to download and align structures.
        :param ensemble_file: (Optional) If specified, each model of the provided PDB path will be one frame of the ensemble.
        '''
        
        self.frames = []
        self.codes=[]
        self.chains = []
        self.workingPath = working_path
        self.code = code

        if ensemble_file != "":
            # if its specified a ensembleFile load it
            frameLines=[]
            for line in FileHandler.getLines(ensemble_file):
                if line.find("MODEL") != -1:
                    dummy, code, chain, num = line.split(" ")
                    self.codes.append(code)
                    self.chains.append(chain)
                    frameLines = []
                elif line.find("ENDMDL") != -1:
                    self.frames.append(structure.Structure(code+chain, from_string = frameLines))
                    frameLines = []
                else:
                    frameLines.append(line)
    
    def toPdbFile(self, path):
        '''
        Save ensemble to PDB file
        
        :param path: filename where to save the ensemble.
        '''
        i=1
        for st in self.frames:
            if i==1: FileHandler.writeLine(path, "")
            line = "MODEL %s %s %s" % (self.codes[i-1], self.chains[i-1],i)
            FileHandler.appendLine(path, line)
            FileHandler.appendLines(path, st.toPdb())
            line = "ENDMDL"
            FileHandler.appendLine(path, line)
            i+=1
        
    def addFrame(self,code,chain,struct):
        '''
        Append a frame to an ensemble.
        
        :param code: Code of the frame to append.
        :param chain: Chain of the frame to append.
        :param struct: Structure object to append.
        '''
        self.codes.append(code)
        self.chains.append(chain)
        self.frames.append(struct)
    
    def getFrame(self,number=0):
        '''
        Get a frame of the ensemble.
        
        :param number: Number of frame to be returned.
        :return: Structure object in the position of the requested frame.
        '''
        if len(self.frames) < number: return None
        else: return self.frames[number]
    
    ##################
    # FoldX Operations
    ##################
    def getTotalEnergy(self, terms=[], consider_waters=False):
        '''
        Computes FoldX total energy for each frame of an ensemble.
        
        :param terms: List of energy terms to return, if not specified, all the energy terms will be returned.
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing the energy terms computed for all structures.
        '''
        
        print( "Computing total energy for ensemble..." )
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=foldxHandler.ENERGY_TERMS)
            terms=foldxHandler.ENERGY_TERMS
        
        for i in trange(len(self.frames), file=sys.stdout):
            retDf.loc[self.codes[i]+"_"+self.chains[i]] = \
                  foldxHandler.getTotalEnergy(self.frames[i], consider_waters)[terms].iloc[0] 
        
        print( "Energy computed." )
        return retDf
    
    def getResiduesEnergy(self, considerWaters=False):
        '''
        Computes FoldX energy for all the residues within each frame of an ensemble object.
        
        :param consider_waters: If True, water molecules within the structure will be taken into account.
        :return: Pandas DataFrame containing total energy of residues as index and structures as columns.
        '''
        
        print( "Computing residue energy for ensemble..." )
        columns = [self.codes[i]+"_"+self.chains[i] for i in range(len(self.codes))]
        retDf = pd.DataFrame(columns=columns)
        
        for i in trange(len(self.frames), file=sys.stdout):
            retDf[self.codes[i]+"_"+self.chains[i]] = \
                  foldxHandler.getResiduesEnergy( self.frames[i], considerWaters)["total"]
        
        print( "Energy computed." )
        return retDf
        
    def mutate(self, mutations, terms=[], generate_mutations_ensemble=False):
        '''
        Generate mutated models with FoldX.
        
        :param mutations: List of mutations in FoldX format to be modeled.
        :param terms: List of energy terms to return, if not specified, all the energy terms will be returned.
        :param generate_mutations_ensemble: If True, an Ensemble object with the mutations and its wild types is returned.
        :return: Tuple with DataFrame of the ddG of the generated models and ensemble of the mutations if this parameter its true.
        '''
        
        print( "Computing mutation(s) %s along ensemble..." % mutations )
        
        # if generate mutated ensemble, generate empty ensembles to load
        if generate_mutations_ensemble:
            trajMut = Ensemble(self.code+"_Mut",self.workingPath+"/%s_Mut/" % mutations.replace(",","_").replace(";","_"))
            trajWT = Ensemble(self.code+"_WT",self.workingPath+"/%s_WT/" % mutations.replace(",","_").replace(";","_"))
        else:
            trajMut = None
            trajWT = None
        
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=foldxHandler.ENERGY_TERMS)
            terms=foldxHandler.ENERGY_TERMS
        
        for i in trange(len(self.frames), file=sys.stdout):
            
            #Add chain to each mutation
            finalMutChain = ""
            for mut in mutations.split(","):
                finalMutChain = mut[0]+self.chains[i]+mut[1:]+","
            finalMutChain = finalMutChain[0:-1]+";"
            
            try:
                ddGsDf, mutModels, wtModels = foldxHandler.getMutants(self.frames[i], finalMutChain, 1)
            except:
                print(self.codes[i]+" failed")
                continue
            
            # if generate mutated ensemble, add models to the ensemples
            if generate_mutations_ensemble:
                trajMut.addFrame(self.codes[i]+"_"+ mutations.replace(",","_"), self.chains[i], mutModels[0] )
                trajWT.addFrame(self.codes[i]+"_"+ mutations.replace(",","_"), self.chains[i], wtModels[0] )
            
            # Append ddGs to return DataFrame
            retDf.loc[self.codes[i]+"_"+self.chains[i]] = ddGsDf.iloc[0]
        
        print( "Mutations computed." )
        return retDf, trajMut, trajWT
    
    def repair(self, fix_residues=[], in_place=True):
        '''
        Minimize and complete sidechains of all structures within an ensemble with FoldX.
        
        :param fix_residues: List of residues in FoldX format to be fixed during repair.
        :param in_place: if true, frames are replaced within the ensemble and None is returned, otherwise a new ensemble is returned.
        :return: Structure object with the repaired structure.
        '''
        print( "Repairing structures along ensemble..." )
        
        ret = None
        for i in trange(len(self.frames), file=sys.stdout):
            repSt = self.frames[i].repair(fix_residues=fix_residues, verbose=False)
            
            if in_place:
                self.frames[i] = repSt
            else:
                ret.addFrame(self.codes[i]+"_Rep", self.chains[i], repSt)
                
        print( "Structures repaired." )
        return ret
    
##########################
#### UNIPROT ENSEMBLE
##########################

class UniprotEnsemble( Ensemble ):
    '''
    Class representing ensembles of structures available for an Uniprot entry.
    
    Available structures will be automatically downloaded, aligned and renumbered (whenever is possible).
    '''
    
    def __init__(self, code, working_path, position=None, just_xray=True, max_resolution=10):
        '''
        Constructor of an uniprot ensemble. 
        
        :param code: Uniprot accession of the ensemble to be created.
        :param working_path: Temporary path to download and align structures.
        :param position: (Optional) if specified, only structures containing the position will be retrieved.
        :param just_xray: (Optional) if True, only x-ray structures will be taken into account.
        :param max_resolution: (Optional) max resolution allowed for retrieved structures.
        '''
        
        super().__init__(code, working_path)
        
        #Add crystals available, considering parameters
        pdbs=[]
        chains = []
        p = Orf(code, replace_existent=True)
        for c in p.getCrystalsWithDetails():
            add = True
            for p in c.property:
                if just_xray and p.type == 'method' and p.value!='X-ray':
                    add=False
                if just_xray and p.type == 'resolution' and float(p.value)>max_resolution:
                    add=False
                if position != None and p.type == 'chains'\
                     and ( int(p.value.split("=")[1].split("-")[0]) > position  \
                        or  int(p.value.split("=")[1].split("-")[1]) < position):
                    add=False 
                if p.type == 'chains':
                    chain=p.value[0]
            
            if add:
                pdbs.append(c.id)
                chains.append(chain)
                
        self.codes=pdbs
        self.chains=chains
        self.rmsd=[]
        
        # Download PDBs Align and make an ensemble
        # TODO: The first on the list will be the reference to align no
        refCode = sorted(self.codes)[0]
        
        print( "Master Structure is %s:" % self.codes[0])
        print( "Aligning to master (total: %s)" % (len(self.codes)) )
        
        ref = structure.Structure(refCode)
        refChain = self.chains[0]
        refSeq = ref.getSequence(refChain)
        refOutPath = self.workingPath+refCode+refChain+".pdb"
        FileHandler.writeLines(refOutPath, ref.toPdb())
        newCodes = []
        newFrames = []
        
        for stPos in trange(1,len(self.codes), file=sys.stdout):
            stCode = self.codes[stPos]
            st = structure.Structure(stCode)
            stChain = self.chains[stPos]
            stSeq = ref.getSequence(self.chains[stPos])
            stOutPath = self.workingPath+stCode+stChain+".pdb"
            FileHandler.writeLines(stOutPath, st.toPdb())
            
            aligned, rms = align(stCode, stOutPath, stSeq, stChain, \
                                 refCode, refOutPath, refSeq, refChain, \
                                 1000)
            
            if aligned:
                newSt = structure.Structure(stCode, stOutPath)
                self.rmsd.append(rms)
                newCodes.append(stCode)
                newFrames.append(newSt)
                
        self.frames = [z for _,z in sorted(zip(self.rmsd, newFrames))]
        self.codes = [z for _,z in sorted(zip(self.rmsd,newCodes))]
        self.rmsd = sorted(self.rmsd)
        
        print("Ensemble built")
        print("Total structures aligned: %s" % len(self.codes) )
        

if __name__ == "__main__":
    
    print("started")
    
    
    outTraj = "/home/lradusky/Downloads/P01112/"
    t = UniprotEnsemble("P01112", outTraj,just_xray=True,max_resolution=1.25)
    t.toPdbFile(outTraj+"traj.pdb")

    outTraj = "/home/lradusky/Downloads/P01112/"
    inTraj = "/home/lradusky/Downloads/P01112/traj.pdb"
    t = Ensemble("P01112", outTraj, inTraj)
     
    retDf, trajMut, trajWT = t.mutate("G13A", generate_mutations_ensemble=True)
     
    print( retDf )
    trajMut.toPdbFile(outTraj+"mutated.pdb")
    trajWT.toPdbFile(outTraj+"WT.pdb")
    
    print( t.getTotalEnergy() )
    t.repair()
    t.toPdbFile(outTraj+"repaired.pdb")
    print( t.getTotalEnergy() )
    print( t.getResiduesEnergy() )
    
    print("done")
    
    