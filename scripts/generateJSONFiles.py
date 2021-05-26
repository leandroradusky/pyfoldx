'''
Created on May 3, 2021
@author: lradusky
'''

import sys
sys.path.append('/home/lradusky/Dropbox/workspacesbg/pyFoldX/')

import json
import datetime
from _collections import defaultdict

from pyfoldx.handlers.urlRetrieveHandler import URLRetrieveHandler
from pyfoldx.handlers.fileHandler import FileHandler
from pyfoldx.handlers.systemHandler import SystemHandler
from pyfoldx.structure.misc import ThreeOne, OneThree
from pyfoldx.structure import Structure

pdbJSON = lambda software, software_version, software_url, \
          pdb_code, chains_annotations, sites_annotations : \
    """ {
          "software_version": "%s",
          "resource_entry_url": "%s",
          "resource_version": "%s",
          "data_resource": "%s",
          "release_date": "%s",
          "pdb_id": "%s",
          "chains": [ %s ] ,
          "sites": [ %s ] ,
          "evidence_code_ontology": [
            {
              "eco_term": "biological system reconstruction (modelling)",
              "eco_code": "ECO_0000088"
            },
            {
              "eco_term": "computational experiment",
              "eco_code": "ECO_0000006"
            }
           ]
        }""" % (software_version, software_url, software_version, software, \
                datetime.date.today().strftime("%d/%m/%Y"), pdb_code, chains_annotations, \
                    ",".join(sites_annotations))

getChainJSON = lambda chain_code, residues_annotations: \
    """ { "chain_label": "%s", %s }""" % (chain_code, residues_annotations)

getResidueJSON = lambda pdb_res_label, aa_type, annotations: \
    """ { "pdb_res_label": "%s",  "aa_type": "%s",  "site_data":  %s  } """ % (pdb_res_label, aa_type, annotations)

generateSiteJSON = lambda site_id, label, source_database, source_accession, description: \
    """
    {
      "site_id": %s,
      "label": "%s",
      "source_database": "%s",
      "source_accession": "%s",
      "additional_site_annotations": { 
                 %s
            }
    }
    """ % (site_id, label, source_database, source_accession, '"description":"%s"' % description if description != "" else "")

def generateFrustrationJsons():

    import pandas as pd

    frustra_path = "/home/lradusky/Dropbox/raduspostdoc/FrustrationMutations/frustraDf.tsv"
    df = pd.read_csv(frustra_path, sep="\t",header=0,index_col=None)

    missense_path = "/home/lradusky/Dropbox/raduspostdoc/pyFoldXpaper/DiseaseMutations/missense3d-benchmarking_PDBeKB.csv"
    df_missense = pd.read_csv(missense_path, sep=",",header=0,index_col=None)
    
    for pdb_code in set(df.pdb):
        site_num = 1
        sites = []
        annotations_chain = []

        for pdb_chain in set(df[df.pdb==pdb_code]["sr_wt_ChainRes"]):

            mutsJson = defaultdict(lambda:[])

            for i, row in df[ (df.pdb==pdb_code) & 
                              (df.sr_wt_ChainRes==pdb_chain) ].iterrows():

                    try:
                        print (pdb_code, pdb_chain, int(row.sr_wt_Res), OneThree[row.sr_wt_AA], OneThree[row.sr_mt_AA], float(row.sr_wt_FrstIndex), float(row.sr_mt_FrstIndex))
                    except:
                        continue
                    
                    try:
                        final_annotation = df_missense[ (df_missense["#PDB"] == pdb_code) & 
                                            (df_missense["#CHAIN"] == pdb_chain) &
                                            (df_missense["#PDBPOS"] == int(row.sr_wt_Res)) &
                                            (df_missense["#RESMUT"] == row.sr_mt_AA) ]["#FINAL_ANNOTATION"].values[0]
                    except: 
                        continue

                    deltaFrustra = round(float(row.sr_wt_FrstIndex) - float(row.sr_mt_FrstIndex),3)
                    mutsJson[(pdb_chain,int(row.sr_wt_Res), row.sr_wt_AA)].append({
                                                "aa_variant":OneThree[row.sr_mt_AA],
                                                "site_id_ref": site_num,
                                                "confidence_classification": "high",
                                                "raw_score": deltaFrustra
                                                })

                    sites.append(generateSiteJSON(site_num, \
                                        final_annotation, \
                                        "missense-3d", \
                                        pdb_code+"_"+pdb_chain+"_"+str(row.sr_wt_Res)+"_"+row.sr_wt_AA+">"+row.sr_mt_AA, \
                                        ""))
                    site_num+=1

            if len(mutsJson) == 0: continue

            # Join the chains jsons
            annotations_res = []
            for (chain, residue, rescode) in mutsJson.keys():

                annotations_mut = []
                for mutJson in mutsJson[(chain, residue, rescode)]:
                    annotations_mut.append(json.dumps(mutJson))
                annotations_res.append(getResidueJSON(residue, OneThree[rescode], json.dumps(mutsJson[(chain,residue,rescode)]))) 

            annotations_chain.append(getChainJSON(pdb_chain,'"residues": ['+",".join(annotations_res)+"]"))

        if len(annotations_chain) == 0: continue

        # The whole PDB json
        json_file = pdbJSON("Frustratometer", "2.0", "http://frustratometer.qb.fcen.uba.ar/", \
                            pdb_code,''+",".join(annotations_chain)+"", sites)
        
        FileHandler.ensureDir("/data/frustra_json/"+pdb_code[1:3]+"/")
        path = "/data/frustra_json/"+pdb_code[1:3]+"/"+pdb_code+".json"
        FileHandler.writeLine(path,  json.dumps(json.loads(json_file), indent=4))


def generateFoldXJsons():
    pass


if __name__ == "__main__":
    
    print( "started" )
    
    generateFrustrationJsons()
    #generateFoldXJsons()

    print( "done" )
