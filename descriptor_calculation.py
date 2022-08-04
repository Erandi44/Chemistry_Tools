'''
Date: 2022-08-03
Author: Erandika Karunaratne
Object: This script calculates following 2D-molecular descriptors for set of molecules provides as an sdf file.

            ‘ExactMolWt’,
            ‘MolLogP’,
            ‘NumHAcceptors’,
            ‘NumHDonors’
            ‘NumAliphaticRings’, 
            ‘NumAliphaticCarbocycles’, 
            ‘NumAliphaticHeterocycles, 
            ‘NumAromaticRings’
            ‘NumAromaticCarbocycles’, 
            ‘NumAromaticHeterocycles’, 
            ‘NumHeteroatoms’, 
            ‘NumRotatableBonds’, 
            ‘NumSaturatedCarbocycles’, 
            ‘NumSaturatedHeterocycles’, 
            ‘NumSaturatedRings’, 
            ‘RingCount’, 
            ‘HeavyAtomCount’, 
            ‘NHOHCount’, 
            ‘NOCount’,
            ‘MaxPartialCharge’, 
            ‘MinPartialCharge’


'''

#import packages
import rdkit
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import display, Image
from rdkit.Chem import Draw
from IPython.core.display import HTML

PandasTools.RenderImagesInAllDataFrames(images=True)




def retrieve_chembl_id(path_input_sdf):
    '''This funtion retrieves PubChem CID no from sdf files, downloaded from PubChem.'''
    cid_no = []
    
    with open(path_input_sdf, "r") as f:
        for line in f:
            #searchphrase = '> <PUBCHEM_COMPOUND_CID>'#for PubChem structures
            searchphrase = '> <chembl_id>'#for chEMBL structures
            
            if searchphrase in line: 
                cid = next(f)
                cid = cid.strip('\n')
                cid_no.append(cid)
            
        print(cid_no)
        return(cid_no)
    
    
    
def calculate_smiles(path_input_sdf):
    '''This funtion retrieves Chembl ID no from sdf files, downloaded from chEMBL.'''
    smiles = []
    suppl = Chem.SDMolSupplier(path_input_sdf, sanitize = False)
    
    
    for mol in suppl:
        try:
            m = Chem.MolToSmiles(mol)
            smiles.append(m)
            print(m)
                        
        except:
            print("An error occured")
            
            
    #chem_df = pd.DataFrame(data=zip(cid_no,smiles),columns=['ChEMBL_ID','SMILES'])
    
    print(smiles)
    return(smiles)


def descriptor_calc(cid_no,smiles,output_path):
    
    '''This function calculates the descriptors and return a csv file with chEMBL_ID/PubChem_ID, 
    SMILES and calculated descriptors'''
    
    #size = (120, 120)
    #mol_image =[]
    
    smiles = smiles
    cid_no = cid_no
    
    #descriptors
    mol_wt = []
    logp =[]
    hba = []
    hbd = []
    aliph_rings =[]
    alph_carbo_cycles = []
    aliph_het_cycles = []
    aromatic_rings = []
    aromatic_carbo_cycles =[]
    aromatic_het_cycles = []
    het_atoms = []
    rotatable_bonds = []
    sat_carbo_cycles = []
    sat_hetero_cycles = []
    sat_rings = []
    ring_count = []
    heavy_atoms = []
    NHOH = []
    NO = []
    max_partial_charge =[]
    min_partial_charge =[]

    
    for m in smiles:
        # create molecular object
        mol = Chem.MolFromSmiles(m)
        
        #Draw 2D structure
        #fig = Draw.MolToMPL(mol, size=size)
        #mol_image.append(fig)
        
        #calculate descriptors
        mol_wt.append(Descriptors.ExactMolWt(mol))
        logp.append(Descriptors.MolLogP(mol))
        hba.append(Descriptors.NumHAcceptors(mol))
        hbd.append(Descriptors.NumHDonors(mol))
        aliph_rings.append(Descriptors.NumAliphaticRings(mol))
        alph_carbo_cycles.append(Descriptors.NumAliphaticCarbocycles(mol))
        aliph_het_cycles.append(Descriptors.NumAliphaticHeterocycles(mol))
        aromatic_rings.append(Descriptors.NumAromaticRings(mol))
        aromatic_carbo_cycles.append(Descriptors.NumAromaticCarbocycles(mol))
        aromatic_het_cycles.append(Descriptors.NumAromaticHeterocycles(mol))
        het_atoms.append(Descriptors.NumHeteroatoms(mol))
        rotatable_bonds.append(Descriptors.NumRotatableBonds(mol))
        sat_carbo_cycles.append(Descriptors.NumSaturatedCarbocycles(mol)) 
        sat_hetero_cycles.append(Descriptors.NumSaturatedHeterocycles(mol))
        sat_rings.append(Descriptors.NumSaturatedRings(mol))  
        ring_count.append(Descriptors.RingCount(mol))
        heavy_atoms.append(Descriptors.HeavyAtomCount(mol))
        NHOH.append(Descriptors.NHOHCount(mol))     
        NO.append(Descriptors.NOCount(mol))  
        max_partial_charge.append(Descriptors.MaxPartialCharge(mol))  
        min_partial_charge.append(Descriptors.MinPartialCharge(mol))
        
        #Draw.MolToMPL(mol, size=(200, 200))
        
    list_collection = list(zip(cid_no, smiles,mol_wt, logp,hba,hbd,aliph_rings,
                               alph_carbo_cycles, aliph_het_cycles,aromatic_rings,
                               aromatic_carbo_cycles, aromatic_het_cycles, het_atoms,
                               rotatable_bonds,sat_carbo_cycles, sat_hetero_cycles,
                               sat_rings, ring_count, heavy_atoms, NHOH, NO,
                               max_partial_charge, min_partial_charge))
                                                     
    desc_cal = pd.DataFrame (list_collection, columns=['chembl_id',
                                                      'SMILES',
                                                      'Mol_Wt',
                                                      'logP',
                                                      'NumHAcceptors',
                                                      'NumHDonors',
                                                      'NumAliphaticRings',
                                                      'NumAliphaticCarbocycles',
                                                      'NumAliphaticHeterocycles',
                                                      'NumAromaticRings',
                                                      'NumAromaticCarbocycles',
                                                      'NumAromaticHeterocycles',
                                                      'NumHeteroatoms',
                                                      'NumRotatableBonds',
                                                      'NumSaturatedCarbocycles',
                                                      'NumSaturatedHeterocycles',
                                                      'NumSaturatedRings',
                                                      'RingCount',
                                                      'HeavyAtomCount',
                                                      'NHOHCount',
                                                      'NOCount',
                                                      'MaxPartialCharge',
                                                      'MinPartialCharge'])
    #print(desc_cal)
    desc_cal.to_csv(output_path+'/'+'descriptor_calculation.csv', index=False)
    
    
def main():
    path_input_sdf = input("Enter Path for input sdf file: ")
    output_path = input("Enter Path for output csv file: ")
    cid_no = retrieve_chembl_id(path_input_sdf)
    smiles = calculate_smiles(path_input_sdf)
    descriptor_calc(cid_no,smiles,output_path)    
    
if __name__ == '__main__':
    main()
        
        

        
