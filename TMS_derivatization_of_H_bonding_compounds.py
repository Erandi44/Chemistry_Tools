'''
Date: 2021/5/20
Objective: This program reads set of compounds from a input sdf file, and derivatize H bond donating groups(i.e. -OH,NH,SH) with TMS (-Si(CH3)3).
Author: Erandika Karunaratne
'''


import rdkit
import os
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit import RDConfig

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcNumHBD




def tms_functionalization(path_input_sdf,cid_list,path_derivatized_sdf):
    #This replace -OH,NH,-SH groups with -Si(CH3)3, compute 2D coordinates, add explicit hydrogens, generate 3D coordinates and optimize coordinates for a given set of molecules
    suppl = Chem.SDMolSupplier(path_input_sdf, sanitize = False)
    
    i=0
    for mol in suppl:
        m = Chem.MolToSmiles(mol)
        start_mol = Chem.MolFromSmiles(m)
        num_hbd = rdkit.Chem.rdMolDescriptors.CalcNumHBD(start_mol)
        mod_mol = [start_mol]
        
        while  num_hbd >0:
            mod_mol = Chem.ReplaceSubstructs(start_mol, Chem.MolFromSmarts('[OX2H]'),Chem.MolFromSmiles(('O[Si](C)(C)C')),replaceAll=True)
            mod_mol = Chem.ReplaceSubstructs(mod_mol[0], Chem.MolFromSmarts('[NX3;H2]'),Chem.MolFromSmiles(('N[Si](C)(C)C')),replaceAll=True)
            mod_mol = Chem.ReplaceSubstructs(mod_mol[0], Chem.MolFromSmarts('[NX3;H]'),Chem.MolFromSmiles(('N[Si](C)(C)C')),replaceAll=True)
            mod_mol = Chem.ReplaceSubstructs(mod_mol[0], Chem.MolFromSmarts('[SX2H]'),Chem.MolFromSmiles(('S[Si](C)(C)C')),replaceAll=True)
            
            Chem.SanitizeMol(mod_mol[0])
            Chem.AssignStereochemistry(mod_mol[0], force=True, cleanIt=True)
            #mod_mol = Chem.MolFromSmiles(mod_mol[0],sanitize=False)
            
            num_hbd = rdkit.Chem.rdMolDescriptors.CalcNumHBD(mod_mol[0])
            #print(num_hbd)
            
            if rdkit.Chem.rdMolDescriptors.CalcNumHBD(mod_mol[0])<=0:
                try:
                    AllChem.Compute2DCoords(mod_mol[0])
                    m3 = Chem.AddHs(mod_mol[0])
                    AllChem.EmbedMolecule(m3,enforceChirality=False,ignoreSmoothingFailures = False, randomSeed=0xf00d)
                    AllChem.MMFFOptimizeMolecule(m3)
                    writer = Chem.SDWriter(path_derivatized_sdf+'/'+cid_list[i]+'.sdf')
                    writer.write(m3)
                    
                    i = i+1   
                    
                except:
                    print("An exception occurred")
                    
                    
                    
def import_sdf(path_input_sdf,path_derivatized_sdf):
    #This extract CID numbers from sdf files and write them into cid_list.text file
    
    dframe =  PandasTools.LoadSDF(path_input_sdf)
    dframe.head()
    cid_list = dframe['PUBCHEM_COMPOUND_CID'].tolist()
    #print(cid_list)
    #len(cid_list)
    
    #write cid numbers of candidates into a text file 
    with open(path_derivatized_sdf+'/cid_list.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % place for place in cid_list)
    
    return(cid_list)


def main():
    path_input_sdf = 'inputpath'  #filepath to input sdf file
    path_derivatized_sdf = 'outputpath'#filepath to a directory of derivatized sdf files (output files)
    cid_list = import_sdf(path_input_sdf,path_derivatized_sdf)
    tms_functionalization(path_input_sdf,cid_list,path_derivatized_sdf)
    
        
if __name__ == "__main__":
    main()
