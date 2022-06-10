'''
Date: 2021/3/15
Objective: This program reads set of compounds from a input sdf file, and generate list of SMILES.
Author: Erandika Karunaratne
'''


import rdkit
import os
from rdkit import Chem



def calculate_SMILES(path_input_sdf,path_output_file,output_filename):
    
    #This converts SDF file to SMILES representation
    suppl = Chem.SDMolSupplier(path_input_sdf, sanitize = False)
    output_file = open (f"{path_output_file}/{output_filename}","w")             
    smiles_list = []
    i=0
    
    for mol in suppl:
        m = Chem.MolToSmiles(mol)
        #print(m)
        i = i+1
        smiles_list.append(m)
                 
                 
    for mol in smiles_list:
        output_file.write(mol + "\n")
    output_file.close()
    
    print(smiles_list)
    #print(len(smiles_list))
    return smiles_list
        
    

def main():
    path_input_sdf = input("Enter file path to the sdf file: ")  #filepath to input sdf file
    path_output_file = input("Enter file path to save the SMILES list: ") #filepath to write SMILES list
    output_filename = input("Enter output file name: ") #file name of SMILES list
    calculate_SMILES(path_input_sdf,path_output_file,output_filename)
    
    
if __name__ == "__main__":
    main()
    
    
    
    