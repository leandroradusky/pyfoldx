'''
Created on Jul 28, 2020
@author: lradusky
'''

from pyfoldx.handlers.FileHandler import FileHandler
from pyfoldx.handlers.URLRetrieveHandler import URLRetrieveHandler
from pyfoldx.handlers.SystemHandler import SystemHandler
from pyfoldx.handlers.XMLToPy import XML2Py
import gzip

SUFFIX = "http://uniprot.org/uniprot"

def obj_dic(d):
    '''
    Convert a dictionary into an object.
    '''
    top = type('new', (object,), d)
    seqs = tuple, list, set, frozenset
    for i, j in d.items():
        if isinstance(j, dict):
            try:
                setattr(top, i, obj_dic(j))
            except:
                pass
        elif isinstance(j, seqs):
            try:
                setattr(top, i, type(j)(obj_dic(sj) if isinstance(sj, dict) else sj for sj in j))
            except:
                pass
        else:
            try:
                setattr(top, i, j)
            except:
                pass
    return top

class Orf(object):
    '''
    Class to handle Uniprot entries, it downloads its xml if not stored in filesystem and convert it to an object.
    '''
    
    def __init__(self,id,replace_existent=False):
        '''
        Constructor of an orf object from its code.
        
        :param id: Uniprot accession (6 letters).
        :param replace_existent: If true, force to re-download the xml file from uniprot.
        '''
        self.id = id
        self.fileName = '/data/uniprot/'+self.id[-3:-1]+'/'+self.id+'.gz' # @TODO: /data has to be parametrizable
        
        # Upload File
        if FileHandler.fileExists(self.fileName) and not replace_existent:
            self.fileLines = FileHandler.getGZLines(self.fileName)
        # If there isn't, download the file
        else:
            self.fileLines = URLRetrieveHandler.RetrieveFileLines('https://www.uniprot.org/uniprot/'+id+'.xml')
            self._saveFile()
        
        xml_lines=""
        with gzip.open(self.fileName, 'r') as f:
            for l in f.readlines():
                xml_lines+= str(l.decode("ascii").strip()).replace(SUFFIX, "")
        
        try:
            uniprot = XML2Py.parse(XML2Py(), xml_lines ) 
        except:
            self.fileLines = URLRetrieveHandler.RetrieveFileLines('http://www.uniprot.org/uniprot/'+id+'.xml')
            self._saveFile()
            xml_lines = map(lambda x: x.replace(SUFFIX, ''), self.fileLines)
            uniprot = XML2Py.parse(XML2Py(), xml_lines )
        
        self.orf = obj_dic(uniprot).uniprot.entry
        
        try:
            accession= self.orf.accession[0].text 
        except: 
            accession= self.orf.accession.text
        
        self.organism = ""
        try:
            for orgname in self.orf.organism.name:
                if orgname.type=='scientific':
                    self.organism = orgname.text
        except:
            if self.orf.organism.name.type=='scientific':
                self.organism = self.orf.organism.name.text
            

        self.accession = accession
        
        self._setName()

    def _saveFile(self):
        FileHandler.writeLines(self.fileName[:-3], self.fileLines)
        SystemHandler.getCommandResult("gzip -f "+self.fileName[:-3])
        
    
    def _setName(self):
        try:
            name= self.orf.protein.submittedName.fullName.text
        except:
            try:
                name= self.orf.protein.recommendedName.fullName.text
            except:
                try:
                    name= self.orf.protein.alternativeName.fullName.text
                except:
                    try:
                        name= self.orf.protein.alternativeName[0].fullName.text
                    except:
                        name = "Putative uncharacterized protein"
                        
        self.name = name 
    
    def getGeneName(self):
        '''
        Extract the gene name of the xml structure.
        
        :return: The gene name of the entry, the primary or the first to appear.
        '''
        try:
            for reference in self.orf.gene.name:
                if reference.type == 'primary':
                    return reference.text
            for reference in self.orf.gene.name:
                return reference.text
        except:
            
            try:
                return self.orf.gene.name.text
            
            except:
                
                return ""

    def getBestCrystal(self):
        '''
        Get the best resolution crystal for this protein.
        
        :return: The best resolution crystal for this protein.
        '''
        crystal = None
        best_res = 100
        # The structures of the protein
        for reference in self.orf.dbReference:
            try:
                if reference.type == 'PDB':
                    for property in reference.property:
                        if property.type == 'resolution':
                            resolution = property.value
                        if property.type == 'chains':
                            chains = property.value
                        if property.type == 'method':
                            method = property.value
                        
                    if method == "X-ray":
                        if crystal == None or best_res > resolution:
                            crystal =reference.id
                            best_res = resolution
            except:
                pass
        return crystal
    
    def getXrayCrystals(self, res_threshold):
        '''
        Get all the x-ray crystals for this protein.
        
        :return: All the x-ray crystals for this protein.
        '''
        crystals = set()
        for reference in self.orf.dbReference:
            try:
                if reference.type == 'PDB':
                    for property in reference.property:
                        if property.type == 'resolution':
                            resolution = property.value
                        if property.type == 'chains':
                            chains = property.value
                        if property.type == 'method':
                            method = property.value
                            
                    if method == "X-ray" and float(resolution) < res_threshold:
                        crystals.add(reference.id)
            except:
                pass
        
        return crystals
        
    
    def getCrystals(self):
        '''
        Get all the crystals for this protein.
        
        :return: All the crystals for this protein.
        '''
        crystals = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'PDB':
                crystals.add(reference.id)
        return crystals
    
    def getCrystalsWithDetails(self):
        '''
        Get all the crystals with its detailed information stored in uniprot.
        
        :return: All the crystals with its detailed information stored in uniprot.
        '''
        crystals = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'PDB':
                crystals.add(reference)
        return crystals
    
    def getPfamFamilies(self):
        '''
        Get all the pfam families for this protein.
        
        :return: All the pfam families for this protein.
        '''
        pfam = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'Pfam':
                pfam.add(reference.id)
        return pfam
    
    def getPfamFamiliesWithDetails(self):
        '''
        Get all the pfam families with its detailed information stored in uniprot.
        
        :return: All the pfam families with its detailed information stored in uniprot.
        '''
        pfam = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'Pfam':
                pfam.add(reference)
        return pfam
    
    def getSNPs(self):
        '''
        Get all the SNPs stored in uniprot.
        
        :return: All the SNPs stored in uniprot as a list of tuples (snp_id,description, original, variation, position).
        '''
        ret = []
        try:
            features = self.orf.feature
        except:
            return ""
        try:
            for f in features:
                pass
        except:
            features =[features]
        
        # The SNPs in the protein
        #try:
        
        for feature in features:
            if feature.type in ["sequence variant", "mutagenesis site"]:
                try:
                    description = feature.description
                except:
                    description = ""
                
                try:
                    original = feature.original.text
                    variation = feature.variation.text
                    position = feature.location.position.position
                except:
                    continue
                try:
                    snp_id = feature.id
                except:
                    snp_id = self.accession+original+variation+position
                    
                ret += [(snp_id,description, original, variation, position)]
    
        return ret
    
if __name__ == "__main__":
    p = Orf('Q9SM56', replace_existent=True)
    for c in p.getCrystalsWithDetails():
        print(c.id)
        for p in c.property:
            print(p.type)
            print(p.value)
        print()
    