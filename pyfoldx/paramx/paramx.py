#!/usr/bin/python
'''
Created on Mar 23, 2020
@author: lradusky
'''

from pyfoldx.handlers.fileHandler import FileHandler

import argparse
from collections import defaultdict
import simplejson

##################################################################################################
#
# CONSTANTS AND GENERAL DEFINITIONS
#
##################################################################################################

# @TODO: make this an option to specify!
ROTABASE_PATH = "rotabase.txt"

#Types
TYPENAMES = defaultdict(str)
TYPENAMES["xyz_sal"] = "water_rotamer"
TYPENAMES["xyz_water"] = "water_rotamer"
TYPENAMES["AAproperties"] = "aaprops"
TYPENAMES["scentropy"] = "tablesidechainentropy"
TYPENAMES["hbondinfo"] = "Hbonds"
TYPENAMES["H_Pos"] = "position_hydrogen"
TYPENAMES["solvenergy"] = "tablesisolvation"

#Fields
FIELDS = defaultdict(str)
FIELDS["water_rotamer"] = ["ion_atom","ion_aa","aa","group","x","y","z","reference_index"]
FIELDS["aaprops"] = ["aa","natural","molweight","extinction","DGmax","movingProtons","virtualProtonsHolder",
                  "numberOfMovingProtons","isAlcohol","protonationStates","misplacedDihed"]
FIELDS["tablesidechainentropy"]=["aa", "entropy", "entropy_abagyan", "radius", "centre_atom", "second_atom"]
FIELDS["Hbonds"] = ["aa","atom","donor","acceptor","number_of_H","number_dummy_H","moving_H",
                    "double_bond","charge","Hname","pKa","dipoled","charged","min_len_hbond",
                    "Hbond_plane_dist","partcov","explicit_solv","tolerance_hbond","hybridization"]
FIELDS["position_hydrogen"]=["aa", "group", "partner1", "pospart1", "partner2", "pospart2", 
                             "hydrogen", "x", "y", "z", "isvirtual", "isexplicit", "protonated", "is_carbonyl"]
FIELDS["tablesisolvation"]=["aa","onelet","atom_type","res_type","atom","volume","Occ","Occmax","vdw",
                            "enerG","level","min_radius","radius","hydrophobic","backbone_atom","polar",
                            "cycle","xTemplate","yTemplate","zTemplate","omega","phi","psi","chi1","chi2",
                            "chi3","chi4","chi5","chi6","chi7","lastDihedral","molecule_type","neigh_atom___16"]

TABLE_ALIAS = defaultdict(str)
TABLE_ALIAS["xyz_sal"] = "saltRotamer"
TABLE_ALIAS["xyz_water"] = "waterRotamer"
TABLE_ALIAS["AAproperties"] = "aminoAcidProperties"
TABLE_ALIAS["scentropy"] = "aminoAcidEntropy"
TABLE_ALIAS["hbondinfo"] = "hydrogenBond"
TABLE_ALIAS["H_Pos"] = "hydrogenPosition"
TABLE_ALIAS["solvenergy"] = "atomSolvation"

ALIAS_TABLE = defaultdict(str)
for k, v in TABLE_ALIAS.items(): ALIAS_TABLE[v] = k 

FIELD_ALIAS = defaultdict(lambda:defaultdict(str))
FIELD_ALIAS["xyz_sal"]["group"] = "atom"
FIELD_ALIAS["xyz_sal"]["ion_atom"] = "ionAtom"
FIELD_ALIAS["xyz_sal"]["ion_aa"] = "ionAminoAcid"

FIELD_ALIAS["xyz_water"]["group"] = "atom"
FIELD_ALIAS["xyz_water"]["ion_atom"] = "ionAtom"
FIELD_ALIAS["xyz_water"]["ion_aa"] = "ionAminoAcid"

FIELD_ALIAS["AAproperties"]["aa"] = "aminoAcid"
FIELD_ALIAS["AAproperties"]["natural"] = "isNatural"
FIELD_ALIAS["AAproperties"]["molweight"] = "molecularWeight"
FIELD_ALIAS["AAproperties"]["extinction"] = "extinctionCoefficient"
FIELD_ALIAS["AAproperties"]["DGmax"] = "maxDeltaG"

FIELD_ALIAS["scentropy"]["entropy_abagyan"] = "moleculeEntropy"
FIELD_ALIAS["scentropy"]["centre_atom"] = "centreAtom"
FIELD_ALIAS["scentropy"]["second_atom"] = "secondAtom"

FIELD_ALIAS["hbondinfo"]["number_of_H"] = "hydrogens"
FIELD_ALIAS["hbondinfo"]["number_dummy_H"] = "dummyHydrogens"
FIELD_ALIAS["hbondinfo"]["moving_H"] = "hydrogensMovility"
FIELD_ALIAS["hbondinfo"]["double_bond"] = "doubleBond"
FIELD_ALIAS["hbondinfo"]["Hname"] = "hydrogenName"
FIELD_ALIAS["hbondinfo"]["dipoled"] = "isDipoled"
FIELD_ALIAS["hbondinfo"]["charged"] = "isCharged"
FIELD_ALIAS["hbondinfo"]["min_len_hbond"] = "minBondDistance"
FIELD_ALIAS["hbondinfo"]["Hbond_plane_dist"] = "planeDistance"
FIELD_ALIAS["hbondinfo"]["partcov"] = "partialCovalentContribution"
FIELD_ALIAS["hbondinfo"]["explicit_solv"] = "explicitSolvation"
FIELD_ALIAS["hbondinfo"]["tolerance_hbond"] = "bondTolerance"

FIELD_ALIAS["position_hydrogen"]["group"] = "atom"
FIELD_ALIAS["position_hydrogen"]["partner1"] = "atomPartner1"
FIELD_ALIAS["position_hydrogen"]["partner2"] = "atomPartner2"
FIELD_ALIAS["position_hydrogen"]["hydrogen"] = "hydrogenName"
FIELD_ALIAS["position_hydrogen"]["isvirtual"] = "isVirtual"
FIELD_ALIAS["position_hydrogen"]["isexplicit"] = "isExplicit"
FIELD_ALIAS["position_hydrogen"]["protonated"] = "isProtonated"
FIELD_ALIAS["position_hydrogen"]["is_carbonyl"] = "isCarbonyl"

FIELD_ALIAS["solvenergy"]["Occ"] = "minOccupancy"
FIELD_ALIAS["solvenergy"]["Occmax"] = "maxOccupancy"
FIELD_ALIAS["solvenergy"]["enerG"] = "solvationEnergy"
FIELD_ALIAS["solvenergy"]["min_radius"] = "vdwInternalRadius"
FIELD_ALIAS["solvenergy"]["radius"] = "vdwClashesRadius"
FIELD_ALIAS["solvenergy"]["hydrophobic"] = "isHydrophobic"
FIELD_ALIAS["solvenergy"]["backbone_atom"] = "isBackbone"
FIELD_ALIAS["solvenergy"]["polar"] = "isPolar"
FIELD_ALIAS["solvenergy"]["cycle"] = "cycleNumber"
FIELD_ALIAS["solvenergy"]["xTemplate"] = "x"
FIELD_ALIAS["solvenergy"]["yTemplate"] = "y"
FIELD_ALIAS["solvenergy"]["zTemplate"] = "z"
FIELD_ALIAS["solvenergy"]["neigh_atom"] = "neighbourAtoms"

##################################################################################################
#
# CONSTANTS AND GENERAL DEFINITIONS
#
##################################################################################################

class JSonParameterLoader( object ):
    
    @staticmethod
    def loadTable(table):
        lines = FileHandler.getLines(ROTABASE_PATH)
        read = False
        
        json_lines=[]
        for line in lines:
            if line == "JSonDataEnd %s" % table or line == "JSonDataEnd %s" % TABLE_ALIAS[table]:
                read = False
                
            if read:
                json_lines.append(line)
            
            if line == "JSonDataStart %s" % table or line == "JSonDataStart %s" % TABLE_ALIAS[table]:
                read = True

        loaded_dict = simplejson.loads("\n".join(json_lines))
        
        return loaded_dict
        
    @staticmethod
    def getAAFieldName(table):
        #if table in ["xyz_water","xyz_sal", "hbondinfo", "H_Pos", "solvenergy", "AAproperties"]:
        return "aa"
    
    @staticmethod
    def getAtomFieldName(table):
        if table in ["xyz_water","xyz_sal","H_Pos"] or ALIAS_TABLE[table] in ["xyz_water","xyz_sal","H_Pos"]:
            return "group"
        elif table in ["hbondinfo", "solvenergy"]:
            return "atom"
        return "atom"
    
    @staticmethod
    def getMoleculesInTable(table):
        loaded_dict = JSonParameterLoader.loadTable(table)
        
        if table not in loaded_dict.keys(): table = TABLE_ALIAS[table]
        
        molecules = set()
        for record in loaded_dict[table]:
            molecules.add(record[JSonParameterLoader.getAAFieldName(table)])
        
        return sorted(list(molecules))
    
    @staticmethod
    def getAtomsInTable(table,molecule):
        loaded_dict = JSonParameterLoader.loadTable(table)
        
        if table not in loaded_dict.keys(): table = TABLE_ALIAS[table]
        
        atoms = set()
        for record in loaded_dict[table]:
            if record[JSonParameterLoader.getAAFieldName(table)] == molecule:
                atoms.add(record[JSonParameterLoader.getAtomFieldName(table)])
        
        return sorted(list(atoms))
    
    @staticmethod
    def getRecordsInTable(table,molecule, atom):
        loaded_dict = JSonParameterLoader.loadTable(TABLE_ALIAS[table])
        
        
        records = []
        for record in loaded_dict[TABLE_ALIAS[table]]:
            
            if record[JSonParameterLoader.getAAFieldName(table)] == molecule and record[JSonParameterLoader.getAtomFieldName(table)] == atom:
                records.append(record)
        
        return records
    
    @staticmethod
    def getAApropertiesRecord(molecule):
        loaded_dict = JSonParameterLoader.loadTable(TABLE_ALIAS["AAproperties"])
        
        for record in loaded_dict[TABLE_ALIAS["AAproperties"]]:
            if record["aa"] == molecule:
                return record
    
    @staticmethod
    def getScentropyRecord(molecule):
        loaded_dict = JSonParameterLoader.loadTable(TABLE_ALIAS["scentropy"])

        for record in loaded_dict[TABLE_ALIAS["scentropy"]]:
            if record["aa"] == molecule:
                return record
    

class JSonMolecule( object ):
    '''
    Molecule parameterized in JSON format.
    '''
    
    def __init__(self, molname=""):
        '''
        Constructor of JSonMolecule
        
        :param molname: Code of the molecule in PDB format (three letters).
        '''
        self.molname = molname
        self.molparams = [ ]
        
        # Create empty tables
        for k in TABLE_ALIAS.keys():
            self.molparams += [defaultdict(lambda:defaultdict([]))]
            self.molparams[-1]["dataType"] = TABLE_ALIAS[k]
            self.molparams[-1][TABLE_ALIAS[k]] = []
        
    def toJson(self, to_file=""):
        '''
        Output the parameterized molecule in JSON format.
        
        :param to_file: (Optional) if specified, the JSON will be written into the path passed to this parameter.
        :return: JSON in string format of to_file parameter not specified.
        '''
        retdict = { "molName": self.molname, "molCode": self.threeLetterCode, "molParams": self.molparams}
        
        if to_file == "":
            return simplejson.dumps(retdict, indent="\t")
        else:
            #return foldxPM_Globals.fileHandler.writeLine(toFile,json.dumps(retdict, indent="\t"))
            return FileHandler.writeLine(to_file,simplejson.dumps(retdict, indent="\t"))
    
    @staticmethod
    def fromJson(file_path):
        '''
        Load parameterized molecule from JSON file.
        
        :param file_path: Path to the file to load the parameters from.
        :return: JSonMolecule object loaded from file.
        '''
        loaded_dict = simplejson.loads("\n".join(FileHandler.getLines(file_path)))
        ret = JSonMolecule(loaded_dict["molName"])
        ret.molparams = loaded_dict["molParams"]
        ret.threeLetterCode = loaded_dict["molCode"]
        
        return ret
    
    def setMolName(self,name):
        self.molname = name
    
    def setThreeLetterCode(self,tlc):
        self.threeLetterCode=tlc
    
    def _transformRecord(self,record,tableName):
        newRecord = defaultdict(str)
        for key in record.keys():
            if key in FIELD_ALIAS[tableName]:
                newRecord[FIELD_ALIAS[tableName][key]] = record[key]
            else:
                newRecord[key] = record[key]
        
        return newRecord
    
    def insertIntoTable(self,tableName,record):
        for table in self.molparams:
            if TABLE_ALIAS[tableName] in table.keys():
                
                # These are uniq
                if tableName in ['AAproperties','scentropy']:
                    table[TABLE_ALIAS[tableName]] = [self._transformRecord(record,tableName)]
                else:
                    table[TABLE_ALIAS[tableName]] += [self._transformRecord(record,tableName)]
    
    def getAtomsInTable(self, table):
        retList = []
        
        for tableName in self.molparams:
            if TABLE_ALIAS[table] in tableName.keys():
                for subtable in tableName[TABLE_ALIAS[table]]:
                    retList += [subtable[JSonParameterLoader.getAtomFieldName(table)]]
        
        return retList

class PdbMolecule( object ):
    
    def __init__(self, path="", fromString=""):
        
        if path != "":
            self.lines = FileHandler.getLines(path)
        else:
            if fromString == "":
                print("Error: Molecule needs to be created from path or string.")
                exit(0)
            else:
                if type(fromString) == type(""):
                    fromString = fromString.split("\n")
                self.lines = fromString
                
        self.atoms = dict()
        self.connections = dict
        self.code = ""
        lineNum = 1
        for line in self.lines:
            if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM':
                atom = dict()
                atom["conect"] = []
                
                # Extract the important fields
                atom["atom_num"] = str(lineNum)
                atom["atom_num_ori"] = line[6:11].strip()
                atom["atom_code"] = line[11:16].strip()
                atom["chain"] = line[21:22].strip()
                atom["seqno"] =line[22:27].strip()
                atom["x"] = float(line[30:38].strip())
                atom["y"] = float(line[38:46].strip())
                atom["z"] = float(line[46:54].strip())
                atom["bfactor"] = float(line[60:66].strip())
                
                #Check is only one one molecule
                if self.code == "" or self.code == line[17:20].strip():
                    self.code = line[17:20].strip()
                else:
                    print( "There is more than one molecule in the provided molecule file." )
                    print( "Exiting..." )
                    exit(0) 
                
                
                self.atoms[atom["atom_code"]] = atom
                lineNum+=1
            
            if line[0:6].strip() == 'CONECT':
                fields = line.split()
                at = self.getAtomNameByNumber(fields[1])
                for field in fields[2:]:
                    self.atoms[at]["conect"].append(self.getAtomNameByNumber(field))
                    
    def getAtomNameByNumber(self, number):
        for atom in self.atoms:
            if self.atoms[atom]["atom_num_ori"] == number:
                return atom
    
##################################################################################################
#
# AUTOPARAMETERIZER MAIN FUNCTIONS
#
##################################################################################################

def checkParameters(args):
    
    if args.molfile == None or not FileHandler.fileExists(args.molfile):
        print( "PDB formatted molecule not provided or not found!" )
        return False
    
    if args.list != None and not FileHandler.fileExists(args.list):
        print( "File with list of atoms not found!" )
        return False
    
    if args.interactive == None:
        print( "Nor interactive mode selected neither list file provided" )
        return False
    
    return True

def loadMoleculeFromPDB(molfile):
    return PdbMolecule(molfile)

def stoleSingleParameter(table, sel_mol, sel_atom, new_mol, new_atom):
    records = JSonParameterLoader.getRecordsInTable(table,sel_mol,sel_atom)
    for record in records:
        record[JSonParameterLoader.getAAFieldName(table)] = new_mol.threeLetterCode
        record[JSonParameterLoader.getAtomFieldName(table)] = new_atom.name
        
        new_mol.insertIntoTable(table, record)

def stoleParameters(molecule, atom):
    
    pass

def checkListFile(mol,listdict):
    
    for atomL in listdict.keys():
        found = False
        for atom in mol.atoms:
            if atomL == mol.atoms[atom]["atom_code"]:
                found = True
                
        if not found:
            print( "Atom %s not found in loaded molecule PDB" % atomL )
            return False
        
        #@TODO: check as well that the atom exists in the molecule of the rotabase
    
    return True

def translateMappings(atom):
    
    if atom == "O_hydroxyl": return ( "OG" , "SER" )
    elif atom =="O_ring": return ( "O4*" , "A" )
    elif atom =="O_minus": return ( "OD1" , "ASP" )
    elif atom =="O_carboxamide": return ( "OD1" , "ASN" )
    elif atom =="N_amino": return ( "NZ" , "LYS" )
    elif atom =="N_guanidino": return ( "NH2" , "ARG" )
    elif atom =="N_imidazol_plus": return ( "ND1" , "HIS" )
    elif atom =="N_imidazol_minus": return ( "NE2" , "HIS" )
    elif atom =="N_pyrazole": return ( "N" , "PRO" )
    elif atom =="N_amide": return ( "ND2" , "ASP" )
    elif atom =="C_ring_not_arom": return ( "CG" , "PRO" )
    elif atom =="C_ring_arom": return ( "CZ" , "PHE" )
    elif atom =="C_single_link": return ( "CG2" , "THR" )
    elif atom =="C_double_link": return ( "CG" , "ARG" )
    elif atom =="C_triple_link": return ( "CG" , "LEU" )
    else:
        print( "The atom with code %s is not recognized. Exiting..." % atom)
        exit(0)

def loadMappings(listfile):
    try:
        retDict = dict()
        for line in FileHandler.getLines(listfile):
            fields = line.split()
            # Assign atom in this molecule (fields[0]) to known atom (fields[1]) of known molecule(fields[2])  
            retDict[fields[0]] = translateMappings(fields[1], fields[2])
    except:
        print( "List with mappings not well formatted, exiting." )
        exit(0)
    
    return retDict
    
def findInListDict(listdict, otheratom):
    for atom in listdict:
        if listdict[atom][0] == otheratom:
            return atom
    return ""
    
def parameterize(mol_string, mappings_dict):
    '''
    Autoparameterize a molecule with its lines and mappings.
    
    :param mol_string: list of PDB lines for the molecule to be parameterized.
    :param mappings_dict: dictionary of mappings to take parameters from. 
    :return: JSonMolecule object parameterized, None if an error arised.
    '''
    
    mol = PdbMolecule(fromString=mol_string)
    JsonMol = JSonMolecule(mol.code)
    JsonMol.setThreeLetterCode(mol.code)
    
    print( "Mappings loaded:" )
    for atom in mappings_dict: 
        if type(mappings_dict[atom]) == type(""):
            mappings_dict[atom] = translateMappings(mappings_dict[atom])
        print( "Atom %s mapped to atom %s" % (atom, mappings_dict[atom] ) )
    
    if not checkListFile(mol,mappings_dict):
        print( "List file error, exiting." )
        return None
    
    # We are ok, we can start to steal :)
    for atom in mappings_dict:
        for table in ["H_Pos", "solvenergy", "hbondinfo", "xyz_water"]:
            # Get the records of the table for the selected molecule and the target atom
            records = JSonParameterLoader.getRecordsInTable(table,mappings_dict[atom][1],mappings_dict[atom][0])
            
            for record in records:
                # Re-define the molecule name in the record
                record[JSonParameterLoader.getAAFieldName(table)] = JsonMol.threeLetterCode
                
                #And here all the magic!
                if table=="solvenergy":
                    record["onelet"]="-"
                    record["res_type"] = "SMALLMOLECULE"
                    record["atom"] = atom
                    record["xTemplate"] = mol.atoms[atom]["x"]
                    record["yTemplate"] = mol.atoms[atom]["y"]
                    record["zTemplate"] = mol.atoms[atom]["z"]
                    
                    record["lastDihedral"]=999
                    
                    record["neigh_atom"] = mol.atoms[atom]["conect"]
                    record["neigh_atom"]+= ["999"] * 16
                    record["neigh_atom"] = record["neigh_atom"][0:16]
                    
                
                if table in ["hbondinfo","xyz_water"]:
                    record["atom"] = atom
                    
                    
                if table=="H_Pos":
                    record["partner1"] = findInListDict(mappings_dict, record["partner1"])
                    record["partner2"] = findInListDict(mappings_dict, record["partner2"])
                    record["group"] = findInListDict(mappings_dict, record["group"])
                    
                    if record["group"] == "" or \
                       record["partner1"] == "" or \
                       record["partner1"] == "":
                        continue
                
                # And insert this record in the corresponding table
                JsonMol.insertIntoTable(table, record)
            
            if len(records) == 0 and table == "xyz_water":
                # DNA ribose is not "waterized", we steal waters from serine
                if mappings_dict[atom][1] in ["A","C","G","T"] and atom[0] == "O":
                    SER_records = JSonParameterLoader.getRecordsInTable(table,"SER","OG")
                    for s_rec in SER_records:
                        record = s_rec
                        record["atom"] = atom
                        record["aa"] = JsonMol.threeLetterCode
                        
                        # And insert this record in the corresponding table
                        JsonMol.insertIntoTable(table, record)
    
    tlc = "SER" #TODO: ask general parameters or keep them fixed? Do they change something in the results?
    try:
        # For the general tables we have to do the same
        record = JSonParameterLoader.getAApropertiesRecord(tlc)
        record["aa"] = JsonMol.threeLetterCode
        JsonMol.insertIntoTable("AAproperties", record)
        record = JSonParameterLoader.getScentropyRecord(tlc)
        record["aa"] = JsonMol.threeLetterCode
        JsonMol.insertIntoTable("scentropy", record)
    except:
        print( "The three letter code that you specified was not found or there was an error." )
        return None
    
    return JsonMol

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--molfile", type=str,
                        help="PDB formatted file with molecule to parameterize")
    parser.add_argument("-i","--interactive", action="store_true", help="Run an atom-by-atom parametrization")
    parser.add_argument("-l","--list", type=str,
                        help="Parameterize molecule with list of atoms to stole from")
    parser.add_argument("-o","--outfile", type=str,
                        help="Output file JSON format")
    args = parser.parse_args()
    
    if checkParameters(args):
        mol = loadMoleculeFromPDB()
        
        print() 
        print( "Read atoms for molecule:" )
        for atom in mol.atoms: 
            print( mol.atoms[atom] )
        print()
        
        mappings = loadMappings(args.list) 
       
        if args.list != None:
            newMol = parameterize(args.molfile, mappings)
        
        print
        if args.outfile == None:
            print( " Produced parameterized file:" )
            print( newMol.toJson() )
        else:
            FileHandler.writeLines(args.outfile, newMol.toJson().split("\n"))
            print( "Parameterized JSON file successfully saved to %s" % args.outfile )
    print() 
    
    print( "Finished Autoparameterizer..." )
    
    