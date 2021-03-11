'''
Created on Nov 8, 2020
@author: lradusky
@summary: miscelaneous definitions and functions
'''

from pyfoldx.handlers.SystemHandler import SystemHandler
import warnings
import sys
import Bio.PDB

# TODO: parameterize this location.
PDB_PATH = '/data/pdb/divided/pdb/'

ThreeOne = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F'
,'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M'
,'ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T'
,'VAL':'V','TRP':'W','TYR':'Y'}

# @summary: align two Structures
def align(code1,file1,seq_str,chain1,code2,file2,use_str,chain2,maxRmsAllowed = 100):
    warnings.filterwarnings("ignore")
    try:
    #if True:
        structure1 = Bio.PDB.PDBParser().get_structure(code1, file1)
        structure2 = Bio.PDB.PDBParser().get_structure(code2, file2)
        
        alt_model = structure1[0]
        ref_chain = structure2[0][chain2]
        alt_chain = structure1[0][chain1]
        
        ref_atoms = []
        alt_atoms = []
        
        lc = lcs(seq_str,use_str)
        lci1 = seq_str.find(lc)
        lci2 = use_str.find(lc)
        
        # Iterate of all chains in the model in order to find all residues
        ref_resl = [j for j in ref_chain][lci2:lci2+len(lc)]
        alt_resl = [j for j in alt_chain][lci1:lci1+len(lc)]
        
        for ref_res, alt_res in zip(ref_resl,alt_resl):
            # Check if residue number ( .get_id() ) is in the list
            try: 
                if ref_res.resname == alt_res.resname:
                    # Append CA atom to list
                    ref_atoms.append(ref_res['CA'])
                    alt_atoms.append(alt_res['CA'])
            except:
                if len(ref_atoms) > len(alt_atoms):
                    ref_atoms.pop(-1)
                pass
        
        #Align these paired atom lists:
        super_imposer = Bio.PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, alt_atoms)
        super_imposer.apply(alt_chain.get_atoms())
        
        #print( "RMS(first model, model %i) = %0.2f" % (alt_model.id, super_imposer.rms) )
        
        if ( super_imposer.rms > maxRmsAllowed) : 
            SystemHandler.executeCommand("rm %s" % file1)
            warnings.filterwarnings("default")
            return False, -1
        
        io=Bio.PDB.PDBIO()
        io.set_structure(alt_model[chain1])
        io.save(file1)
    
    except Exception as e:
        #print("Structures %s and %s were not aligned" % (code1,code2))
        SystemHandler.executeCommand("rm %s" % file1)
        warnings.filterwarnings("default")
        return False, -1
    
    warnings.filterwarnings("default")
    return True, super_imposer.rms

# @summary: largest common string,for alignment
def lcs(S,T):
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = ""
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    longest = c
                    lcs_set = S[i-c+1:i+1]
                

    return lcs_set


import os
import threading
import time


class OutputGrabber(object):
    """
    Class used to grab standard output or another stream.
    """
    escape_char = "\b"

    def __init__(self, stream=None, threaded=False):
        self.origstream = stream
        self.threaded = threaded
        if self.origstream is None:
            self.origstream = sys.stdout
        self.origstreamfd = self.origstream.fileno()
        self.capturedtext = ""
        # Create a pipe so the stream can be captured:
        self.pipe_out, self.pipe_in = os.pipe()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, type, value, traceback):
        self.stop()

    def start(self):
        """
        Start capturing the stream data.
        """
        self.capturedtext = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.origstreamfd)
        # Replace the original stream with our write pipe:
        os.dup2(self.pipe_in, self.origstreamfd)
        if self.threaded:
            # Start thread that will read the stream:
            self.workerThread = threading.Thread(target=self.readOutput)
            self.workerThread.start()
            # Make sure that the thread is running and os.read() has executed:
            time.sleep(0.01)

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # Print the escape character to make the readOutput method stop:
        self.origstream.write(self.escape_char)
        # Flush the stream to make sure all our data goes in before
        # the escape character:
        self.origstream.flush()
        if self.threaded:
            # wait until the thread finishes so we are sure that
            # we have until the last character:
            self.workerThread.join()
        else:
            self.readOutput()
        # Close the pipe:
        os.close(self.pipe_in)
        os.close(self.pipe_out)
        # Restore the original stream:
        os.dup2(self.streamfd, self.origstreamfd)
        # Close the duplicate stream:
        os.close(self.streamfd)

    def readOutput(self):
        """
        Read the stream data (one byte at a time)
        and save the text in `capturedtext`.
        """
        while True:
            char = os.read(self.pipe_out,1).decode(self.origstream.encoding)
            if not char or self.escape_char in char:
                break
            self.capturedtext += char
    
def in_notebook():
    """
    Returns True if the module is running in IPython kernel,
    False if in IPython shell or other Python shell.
    """
    return 'ipykernel' in sys.modules
