--------------------------------------------------------------------------------------------------------------
'''
Date: 2022/06/30

Objective: This program reads a set of compounds from an SDF file and transforms functional groups in 
a molecule based on the requirement of the user.

Author: Erandika Karunaratne
'''
--------------------------------------------------------------------------------------------------------------

import rdkit
import os
import pandas as pd

from rdkit.Chem import PandasTools
from rdkit import RDConfig
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw

from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions

DrawingOptions.bondLineWidth=1.8
DrawingOptions.atomLabelFontSize=14
DrawingOptions.includeAtomNumbers=False

from __future__ import print_function
# import rdkit components

IPythonConsole.ipython_useSVG=True

# for flattening tuples and lists
from itertools import chain



# This program can perform following functional group transformations.
# Improvements needed: 
#    1. Include reactions related to ethers and acid anhydrides.
#    2. Alkene addition postion based on  bulkiness of R group.
#    3. Include all products from products tuple and remove duplicates. 


reactions = f"""Alkyl Halides
                1. alkyl halide(R-X) -> alkyl cyanaide(R-CN)
                2. alkyl halide(R-X) -> alkene (C=C)
                 
                Alkene
                3. alkene(C=C) -> alcohol(R-OH)
                4. alkene(C=C) -> epoxide
                5. alkene(C=C) -> alkyl halide(R-X)

                Alcohol
                6. alcohol(R-OH) -> alkene
                7. primary alcohol(R-OH) -> aldehyde(R-COH)
                8. primary alcohol(R-OH) -> carboxylic acid(R-COOH)
                9. secondary alcohol(R-OH) -> ketone(R-CO-R)

                Aldehyde
                10. aldehyde -> alcohol(R-OH)
                11. aldehyde -> carboxylic acid(R-COOH)

                Ketone
                12. ketone(R-CO-R) -> secondary alcohol(R-OH)

                Carboxylic Acid
                13. carboxylic acid -> primary alcohol
                14. carboxylic acid -> acyl chloride
                15. carboxylic acid(R-COOH) -> ester(R-CO-O-CH3)

                Amide
                16. amide(R-CO-NH2) -> primary amine (R-NH2)
                17. amide(R-CO-NH2) -> alkyl cyanide (R-CN)
                18. amide(R-CO-NH2) -> carboxylic acid (R-COOH)

                Acyl Chloride
                19. acyl chloride(R-CO-Cl) -> carboxylic acid(R-COOH)
                20. acyl chloride(R-CO-Cl) -> amide(R-CO-NH2)
                21. acyl chloride(R-CO-Cl) -> ester(R-CO-O-CH3)

                Ester
                22. ester(R-CO-O-R) -> carboxylic acid (R-COOH)

                Alkyl Cyanide
                23. alkyl cyanide (R-CN) -> primary amine (R-NH2)
                24. alkyl cyanide (R-CN) -> amide(R-CO-NH2)"""



# This list represnt the functional group that needs to be transformed.
fl_gps = ['[C:1][F,Cl,Br,I:2]',
'[C!H0:1][C:2][Cl:3]',
'[C:1]=[C:2]',
'[C:1]=[C:2]',
'[C:1]=[C:2]',
'[C:1]([C:2][OH:3])',
'[C:1][OH:2]',
'[C:1][OH:2]',
'[C:1]-[C:2](-[C:4])-[OH:3]',
'[CX3H1:1]=[O:2]',
'[CX3H1:1]=[O:2]',
'[C:1]-[C:2](-[C:4])=[O:3]',
'[CX3:1](=O)[OX2H1:3]',
'[CX3:1](=O)[OX2H1:3]',
'[CX3:1](=O)[OX2H1:3]',
'[CX3:1](=O)[NX3;H2:2]',
'[CX3:1](=O)[NX3;H2:2]',
'[CX3:1](=O)[NX3;H2:2]',
'[CX3:1](=O)[Cl:3]',
'[CX3:1](=O)[Cl:3]',
'[CX3:1](=O)[Cl:3]',
'[CX3:1](=O)[OX2:3][CX4]',
'[C:1]#[N:2]',
'[C:1]#[N:2]']

# This list reperesent Reaction SMARTS for above chemical transformations.
rxns = ['[C:1][F,Cl,Br,I:2]>>[C:1][C:2]#[N:3]',
'[C!H0:1][C:2][Cl:3]>>[C:1]=[C:2].[Cl:3]',
'[C:1]=[C:2]>>[C:1]([C:2][OH:3])',
'[C:1]=[C:2]>>[OX2r3:3]1[CX4r3:1][CX4r3:2]1',
'[C:1]=[C:2]>>[C:1]([C:2][Cl:3])',
'[C:1]([C:2][OH:3])>>[C:1]=[C:2]',
'[C:1][OH:2]>>[CX3H1:1]=[O:2]',
'[C:1][OH:2]>>[CX3:1](=O)[OX2H1:3]',
'[C:1]-[C:2](-[C:4])-[OH:3]>>[C:1]-[C:2](-[C:4])=[O:3]',
'[CX3H1:1]=[O:2]>>[C:1]([C:2][OH:3])',
'[CX3H1:1]=[O:2]>>[CX3:1](=O)[OX2H1:3]',
'[C:1]-[C:2](-[C:4])=[O:3]>>[C:1]-[C:2](-[C:4])-[OH:3]',
'[CX3:1](=O)[OX2H1:3]>>[C:1]([C:2][OH:3])',
'[CX3:1](=O)[OX2H1:3]>>[CX3:1](=O)[Cl:3]',
'[CX3:1](=O)[OX2H1:3]>>[CX3:1](=O)[OX2:3][CH3:4]',
'[CX3:1](=O)[NX3;H2:2]>>[CX4;H2:1][NX3;H2]',
'[CX3:1](=O)[NX3;H2:2]>>[C:1]#[N:2]',
'[CX3:1](=O)[NX3;H2:2]>>[CX3:1](=O)[OX2H1:2]',
'[CX3:1](=O)[Cl:3]>>[CX3:1](=O)[OX2H1:3]',
'[CX3:1](=O)[Cl:3]>>[CX3:1](=O)[NX3;H2:2]',
'[CX3:1](=O)[Cl:3]>>[CX3:1](=O)[OX2:3][CH3:4]',
'[CX3:1](=O)[OX2:3][CX4]>>[CX3:1](=O)[OX2H1:3]',
'[C:1]#[N:2]>>[CX4;H2:1][NX3;H2:2]',
'[C:1]#[N:2]>>[CX3:1](=O)[NX3;H2:2]'    
]


def user_input(rxns):
    ''' Process user requirement'''
    print("You can perform following functional group trasformations.")
    print("\n")
    print(reactions)    
    user = int(input("Please select a number from above list (between 1-29):"))
    
    if (1 > user >24 ): 
        print("Invalid input")
        user = input("Please select a number from above list (between 1-29):")
    
    user = user-1
    return user



def retrieve_cid_identifiers(input_path):
    '''This funtion retrieves PubChem CID no from sdf files, downloaded from PubChem.'''
    cid_no = []
    
    with open(input_path, "r") as f:
        for line in f:
            searchphrase1 = '> <PUBCHEM_COMPOUND_CID>'
            
            if searchphrase1 in line: 
                cid = next(f)
                cid = cid.strip('\n')
                cid_no.append(cid)
            
        print(cid_no)
        return(cid_no)
    

def molecule_list(input_path):
    '''This function reads sdf file and extract molecules into a list as objects'''
    suppl = Chem.SDMolSupplier(input_path)
    mol_list=[mol for mol in suppl]
    print(mol_list)    
    return mol_list


def chem_transformation(output_path, mol_list,rxns,cid_no,user,fl_gps):
    '''This function performs the chemical transformation and save transformed/not-transformed molecules in a user defined directory.'''
    
    
    # Retrieve functional gp to be transformed
    fl_gp = fl_gps[user]
    patt = Chem.MolFromSmarts(fl_gp) # Assigning SMARTS representation of functional needs to replaced into patt
    
    rxn = AllChem.ReactionFromSmarts(rxns[user])# Extract the chemical reaction
        
    i=0
        
    for mol in mol_list:
        mol.UpdatePropertyCache()
        presence = int(mol.HasSubstructMatch(patt)) # Checking the presence of functional group
        
        # Perform chemical transformation if interested functional group present in the molecule
        # Then sanitize,add H's, generate 3D coordinates, optimize and save as sdf file
        if presence ==1:
            product = rxn.RunReactants((mol, ))
            
            try:
                Chem.SanitizeMol(product[0][0])
                Chem.AssignStereochemistry(product[0][0], force=True, cleanIt=True)
                AllChem.Compute2DCoords(product[0][0])
                mod_mol_with_H = Chem.AddHs(product[0][0]) # add H's
                AllChem.EmbedMolecule(mod_mol_with_H,randomSeed=0xf00d)# convert to 3D coordinates
                AllChem.MMFFOptimizeMolecule(mod_mol_with_H) # optimization
                #print(Chem.MolToMolBlock(mod_mol_with_H))
                writer = Chem.SDWriter(output_path+cid_no[i]+'.sdf')
                writer.write(mod_mol_with_H)
                #print(Chem.MolToMolBlock(mod_mol_with_H))
            except:
                print("An exception occured")
        
        # This will sanitize,add H's, generate 3D coordinates, optimize and save the molecule as a sdf file
        # if it doesn't contain the funtional group of interest.
        else:
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol,force=True, cleanIt=True)
            AllChem.Compute2DCoords(mol)
            mod_mol_with_H = Chem.AddHs(mol) # add H's
            AllChem.EmbedMolecule(mod_mol_with_H,randomSeed=0xf00d)# conver to 3D coordinates
            AllChem.MMFFOptimizeMolecule(mod_mol_with_H) # optimization
            writer = Chem.SDWriter(output_path+cid_no[i]+'.sdf')
            writer.write(mod_mol_with_H)
            #print(Chem.MolToMolBlock(mol))   
        i=i+1
    
        
    
    
def main():
    input_path = input("Enter Path for input sdf file: ")
    output_path = input("Enter Path for output sdf file: ")
    user = user_input(rxns)
    mol_list = molecule_list(input_path)
    cid_no = retrieve_cid_identifiers(input_path)
    chem_transformation(output_path, mol_list,rxns,cid_no,user,fl_gps)
    
    
    
if __name__ == "__main__":
    main()