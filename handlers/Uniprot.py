'''
Created on Jul 28, 2020
@author: lradusky
'''

from src.handlers.FileHandler import FileHandler
from src.handlers.URLRetrieveHandler import URLRetrieveHandler
from src.handlers.SystemHandler import SystemHandler
from src.handlers.XMLToPy import XML2Py
import time
import gzip

SUFFIX = "http://uniprot.org/uniprot"
def obj_dic(d):
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
    
    def __init__(self,id,replaceExistent=False):
        self.id = id
        self.fileName = '/data/uniprot/'+self.id[-3:-1]+'/'+self.id+'.gz'
        
        # Upload File
        if FileHandler.fileExists(self.fileName) and not replaceExistent:
            self.fileLines = FileHandler.getGZLines(self.fileName)
        # If there isn't, download the file
        else:
            self.fileLines = URLRetrieveHandler.RetrieveFileLines('https://www.uniprot.org/uniprot/'+id+'.xml')
            self.saveFile()
        
        xml_lines=""
        with gzip.open(self.fileName, 'r') as f:
            for l in f.readlines():
                xml_lines+= str(l.decode("ascii").strip()).replace(SUFFIX, "")

        try:
            uniprot = XML2Py.parse(XML2Py(), xml_lines ) 
        except:
            self.fileLines = URLRetrieveHandler.RetrieveFileLines('http://www.uniprot.org/uniprot/'+id+'.xml')
            self.saveFile()
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
        
        self.setName()
    
    def getECNumbers(self):
        ret = []
        for reference in self.orf.dbReference:
            if reference.type == 'EC':
                ret.append(reference.id)
        return ret 
    
    def getOntologies(self):
        ret = []
        for reference in self.orf.dbReference:
            if reference.type == 'GO':
                ret.append(reference.id)
        return ret 
        
    def getTaxonomyLineage(self):
        ret = []
        for taxon in self.orf.organism.lineage.taxon:
            ret += [taxon.text]
        return ret
    
    def getDNASeq(self):
        # @todo: actually doesnt work!!
        ret = []
        i=0
        for reference in self.orf.dbReference:
            if reference.type == 'EMBL':
                i+=1
                if i == 2:
                    ret+= [URLRetrieveHandler.RetrieveFileLines("http://www.ebi.ac.uk/ena/data/view/"+reference.id+"&display=fasta")]
                    return ret
        return ret

    def saveFile(self):
        FileHandler.writeLines(self.fileName[:-3], self.fileLines)
        #f_in = open(self.fileName[:-3], 'rb')
        #f_out = gzip.open(self.fileName, 'wb')
        #f_out.writelines(f_in)
        #f_out.close()
        #f_in.close()
        SystemHandler.getCommandResult("gzip -f "+self.fileName[:-3])
        
    
    def setName(self):
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
    
    def getTaxonomy(self):
        return self.orf.organism.dbReference.id

    def getGeneNameOrderedLocus(self):
        try:
                
            # The structures of the protein
            for reference in self.orf.gene.name:
                if reference.type == 'ordered locus':
                    return reference.text
        except:
            return ""
    
    def getKeggCode(self):
        for reference in self.orf.dbReference:
            try:
                if reference.type == 'KEGG':
                    return reference.id
            except:
                pass
        return ""
    
    
    def getGeneName(self):
        try:
                
            # The structures of the protein
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
        crystals = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'PDB':
                crystals.add(reference.id)
        return crystals
    
    def getCrystalsWithDetails(self):
        crystals = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'PDB':
                crystals.add(reference)
        return crystals
    
    
    def getPfamFamilies(self):
        pfam = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'Pfam':
                pfam.add(reference.id)
        return pfam
    
    def getPfamFamiliesWithDetails(self):
        pfam = set()
        # The structures of the protein
        for reference in self.orf.dbReference:
            if reference.type == 'Pfam':
                pfam.add(reference)
        return pfam
    
    
    def getGlycosylationSites(self):
        try:
            features = self.orf.feature
        except:
            return ""
        try:
            for f in features:
                pass
        except:
            features =[features]
        
        ret = []
        
        for feature in features:
            if feature.type in ["glycosylation site"]:
                ret+=[(feature.location.position.position, feature.description)]
        
        return ret
    
    def getSNPs(self):
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
    p = Orf('Q9SM56', replaceExistent=True)
    for c in p.getCrystalsWithDetails():
        print(c.id)
        for p in c.property:
            print(p.type)
            print(p.value)
        print()
    