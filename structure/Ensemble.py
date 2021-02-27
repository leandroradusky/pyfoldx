'''
Created on Nov 2, 2020
@author: lradusky
'''

from src.handlers.FileHandler import FileHandler
import src.structure.Structure as Structure
from src.structure.misc import align
from src.handlers.Uniprot import Orf
from src.pyFoldX import pyFoldX

from glob import glob
import pandas as pd
from tqdm import tqdm,trange
import time

##########################
#### ENSEMBLE
##########################

class Ensemble(object):

    def __init__(self, code, workingPath,ensembleFile=""):
        
        self.frames = []
        self.codes=[]
        self.rmsd=[]
        self.chains = []
        self.workingPath = workingPath
        self.code = code

        if ensembleFile != "":
            # if its specified a ensembleFile load it
            frameLines=[]
            for line in FileHandler.getLines(ensembleFile):
                if line.find("MODEL") != -1:
                    dummy, code, chain, num = line.split(" ")
                    self.codes.append(code)
                    self.chains.append(chain)
                    frameLines = []
                elif line.find("ENDMDL") != -1:
                    self.frames.append(Structure.Structure(code+chain, fromString = frameLines))
                    frameLines = []
                else:
                    frameLines.append(line)
    
    def saveToPDB(self, path):
        i=1
        for st in self.frames:
            if i==1: FileHandler.writeLine(path, "")
            line = "MODEL %s %s %s" % (self.codes[i-1], self.chains[i-1],i)
            FileHandler.appendLine(path, line)
            FileHandler.appendLine(path, st.toPdb())
            line = "ENDMDL"
            FileHandler.appendLine(path, line)
            i+=1
        
    def addFrame(self,code,chain,struct):
        self.codes.append(code)
        self.chains.append(chain)
        self.frames.append(struct)
    
    def getFrame(self,number=0):
        if len(self.frames) >= number: return None
        else: return self.frames[number]
    
    ##################
    # FoldX Operations
    ##################
    def getTotalEnergy(self, terms=[],considerWaters=False):
        
        print( "Computing total energy for ensemble..." )
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=pyFoldX.ENERGY_TERMS)
            terms=pyFoldX.ENERGY_TERMS
        
        for i in trange(len(self.frames)):
            retDf.loc[self.codes[i]+"_"+self.chains[i]] = \
                  pyFoldX.getTotalEnergy(self.codes[i]+"_"+self.chains[i], self.frames[i].toPdb(), considerWaters)[terms].loc[self.codes[i]+"_"+self.chains[i]] 

        
        print( "Energy computed." )
        return retDf
    
    def getResiduesEnergy(self, considerWaters=False):
        
        print( "Computing residue energy for ensemble..." )
        columns = [self.codes[i]+"_"+self.chains[i] for i in range(len(self.codes))]
        retDf = pd.DataFrame(columns=columns)
        
        for i in trange(len(self.frames)):
            retDf[self.codes[i]+"_"+self.chains[i]] = \
                  pyFoldX.getResiduesEnergy( self.frames[i].toPdb(), considerWaters)["total"]
        
        retDf.index = retDf.index.map(lambda x: (("000000"+x[3:])[-5:]) +"_"+ x[0:3])
        retDf.sort_index(inplace=True)
        
        print( "Energy computed." )
        return retDf
        
    def mutate(self, mutations, terms=[], generateMutationsEnsemble=False):
        '''
        :todo: implement several runs and keep the best
        '''
        print( "Computing mutation(s) %s along ensemble..." % mutations )
        
        # if generate mutated ensemble, generate empty ensembles to load
        if generateMutationsEnsemble:
            trajMut = Ensemble(self.code+"_Mut",self.workingPath+"/%s_Mut/" % mutations.replace(",","_").replace(";","_"))
            trajWT = Ensemble(self.code+"_WT",self.workingPath+"/%s_WT/" % mutations.replace(",","_").replace(";","_"))
        else:
            trajMut = None
            trajWT = None
        
        if terms != []:
            retDf = pd.DataFrame(columns=terms)
        else:
            retDf = pd.DataFrame(columns=pyFoldX.ENERGY_TERMS)
            terms=pyFoldX.ENERGY_TERMS
        
        for i in trange(len(self.frames)):
            
            #Add chain to each mutation
            finalMutChain = ""
            for mut in mutations.split(","):
                finalMutChain = mut[0]+self.chains[i]+mut[1:]+","
            finalMutChain = finalMutChain[0:-1]+";"
            
            try:
                ddGsDf, mutModels = pyFoldX.mutate(self.frames[i].toPdb(), finalMutChain, 1)
            except:
                print(self.codes[i]+" failed")
                continue
            
            # if generate mutated ensemble, add models to the ensemples
            if generateMutationsEnsemble:
                trajMut.addFrame(self.codes[i]+"_"+ mutations.replace(",","_"),
                                 self.chains[i], Structure.Structure(self.codes[i]+"_"+ mutations.replace(",","_"), fromString=mutModels[0][0].split("\n")))
                trajWT.addFrame(self.codes[i]+"_"+ mutations.replace(",","_"),
                                 self.chains[i], Structure.Structure(self.codes[i]+"_"+ mutations.replace(",","_"), fromString=mutModels[0][1].split("\n")))
            
            retDf.loc[self.codes[i]+"_"+self.chains[i]] = ddGsDf.loc[0]
            
            
            # Append ddGs to return DataFrame
        
        print( "Mutations computed." )
        return retDf, trajMut, trajWT

    
##########################
#### UNIPROT ENSEMBLE
##########################

class uniprotEnsemble( Ensemble ):
    
    def __init__(self, code, workingPath, position=None, justXray=True, maxResolution=10):
        super().__init__(code, workingPath)
        
        #Add crystals available, considering parameters
        pdbs=[]
        chains = []
        p = Orf(code, replaceExistent=True)
        for c in p.getCrystalsWithDetails():
            add = True
            for p in c.property:
                if justXray and p.type == 'method' and p.value!='X-ray':
                    add=False
                if justXray and p.type == 'resolution' and float(p.value)>maxResolution:
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
        
        # Download PDBs Align and make an ensemble
        # TODO: The first on the list will be the reference to align no
        refCode = self.codes[0]
        
        print( "Master Structure is %s:" % self.codes[0])
        print( "Aligning to master (total: %s)" % (len(self.codes)) )
        
        ref = Structure.Structure(refCode)
        refChain = self.chains[0]
        refSeq = ref.getSequence(refChain)
        refOutPath = self.workingPath+refCode+refChain+".pdb"
        FileHandler.writeLine(refOutPath, ref.toPdb())
        newCodes = []
        newFrames = []
        
        for stPos in range(1,len(self.codes)):
            stCode = self.codes[stPos]
            st = Structure.Structure(stCode)
            stChain = self.chains[stPos]
            stSeq = ref.getSequence(self.chains[stPos])
            stOutPath = self.workingPath+stCode+stChain+".pdb"
            FileHandler.writeLine(stOutPath, st.toPdb())
            
            aligned, rms = align(stCode, stOutPath, stSeq, stChain, \
                                 refCode, refOutPath, refSeq, refChain, \
                                 1000)
            
            if aligned:
                newSt = Structure.Structure(stCode, stOutPath)
                self.rmsd.append(rms)
                newCodes.append(stCode)
                newFrames.append(newSt)
                
            print("%s..." % (stPos), end='', flush=True)
            if stPos%5==0:print("of %s structures processed" % len(self.codes))
            
        self.frames = [z for _,z in sorted(zip(self.rmsd, newFrames))]
        self.codes = [z for _,z in sorted(zip(self.rmsd,newCodes))]
        self.rmsd = sorted(self.rmsd)
        
        print("Ensemble built")
        print("Total structures aligned: %s" % len(self.codes) )
        

##########################
#### PDBFLEX ENSEMBLE
##########################

class pdbFlexEnsemble( Ensemble ):
    
    def __init__(self, code, workingPath, path=""):
        super().__init__(code, workingPath)
        
        self.path = path
        self._processConfigurations()
        self._buildEnsemble()
    
    def _buildEnsemble(self):
        
        newCodes = []
        newChains= []
        newFrames = []
        
        print( "MasterStructure is %s:" % self.code)
        print( "Aligning to master (total: %s)" % (len(self.frames)) )
        
        pTempl =self.workingPath+self.code+".pdb"
        stTempl = Structure.Structure(self.code, pTempl)
        
        j=0
        for c, frame in zip(self.codes, self.frames):
            j+=1
            
            p = self.workingPath+c+".pdb"
            
            seq1 = frame.getSequence(c[-1])
            seq2 = stTempl.getSequence(self.code[-1])
            
            aligned, rms = align(c, p, seq1, c[-1], self.code, pTempl, seq2, self.code[-1], 1000)
            if aligned:
                newSt = Structure.Structure(c, p)
                self.rmsd.append(rms)
                newCodes.append(c)
                newFrames.append(newSt)
                newChains.append(c[-1])
                
            print("%s..." % (j), end='', flush=True)
            if j%10==0:print()
            
        self.frames = [z for _,z in sorted(zip(self.rmsd, newFrames))]
        self.codes = [z for _,z in sorted(zip(self.rmsd,newCodes))]
        self.chains = [z for _,z in sorted(zip(self.rmsd,newChains))]
        self.rmsd = sorted(self.rmsd)
        
        print("Ensemble built")
    
    def _processConfigurations(self):
        pdbs = glob("%s/*.pdb" % self.path)
        print( "Processing configurations (total: %s)" % (len(pdbs)) )
        print( "Processed: " )
        j=0
        for pdb in pdbs:
            j+=1
            
            stTempl = Structure.Structure(pdb, pdb)
            newPdb = pdb.split("/")[-1].split(".")[0][0:4]
            chain = pdb.split("/")[-1].split(".")[0][4]
            try: st = Structure.Structure(newPdb)
            except: continue
            
            outPath = self.workingPath+newPdb+chain+".pdb"
            FileHandler.writeLine(outPath, st.toPdb())
            
            seq1 = st.getSequence(chain)
            seq2 = stTempl.getSequence(chain)
            
            if align(newPdb, outPath, seq1, chain, pdb, pdb, seq2, chain)[0]:
                newSt = Structure.Structure(newPdb+chain, outPath)
                self.frames.append(newSt)
                self.codes.append(newPdb+chain)
            
            print("%s..." % (j), end='', flush=True)
            if j%10==0:print()
            
        print("Configurations processed")

if __name__ == "__main__":
    
    print("started")
    '''
    code = "4k81B"
    sampleTraj = "/home/lradusky/Downloads/%s/" % code
    outTraj = "/home/lradusky/Downloads/outputTraj/"
    
    t = pdbFlexEnsemble(code, outTraj, sampleTraj)
    print( t.rmsd )
    t.saveToPDB(outTraj+"traj.pdb")
    t.mutate(t.getModel().data["B"][22], "V")
    '''
    outTraj = "/home/lradusky/Downloads/P01112/"
    t = uniprotEnsemble("P01112", outTraj,justXray=True,maxResolution=1.6)
    t.saveToPDB(outTraj+"traj.pdb")
    
    outTraj = "/home/lradusky/Downloads/P01112/"
    inTraj = "/home/lradusky/Downloads/P01112/traj.pdb"
    t = Ensemble("P01112", outTraj, inTraj)
    
    retDf, trajMut, trajWT = t.mutate("G13A", generateMutationsEnsemble=True)
    
    trajMut.saveToPDB(outTraj+"mutated.pdb")
    trajWT.saveToPDB(outTraj+"WT.pdb")
    
    t.saveToPDB(outTraj+"traj2.pdb")
    totalEnergyDf = t.getTotalEnergy()
    residuesEnergyDf = t.getResiduesEnergy()
    print( residuesEnergyDf )
    
    print("done")
    
    