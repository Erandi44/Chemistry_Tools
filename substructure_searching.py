'''
Date: 2022/6/14
Objective: This program reads set of compounds from a input sdf file, check for the presence of different functional groups/sub-structures and writes the presence(1) or absence(0) of each functional group into a csv file.
Author : Erandika Karunaratne
'''


import rdkit
import os
import csv
from rdkit import Chem
from rdkit.Chem import AllChem



def sub_structures():

    '''
    Functional groups/substructures that can be searched using this script. Edit this SMARTS list as necessary

    Methyl group [CX4H3]
    Methylene group [CX4H2], [CX4H]
    Alkene	[#6]=[#6]
    Alkyne	[CX2]#[CX2]
    Allene	[CX3]=[CX2]=[CX3]
    chloride	[ClX1]
    fluoride	[FX1]
    bromide	[BrX1]
    iodide	[IX1]
    Phenol	[OX2H][c]
    Alcohol [OX2H]
    Thiol	[SX2H]
    ether	[#6]-[#8]-[#6]
    thioeher	[SX2]
    primary amine	[NX3;H2]
    secondary amine	[NX3;H]
    Disulfide	[SX2D2][SX2D2]
    Aldehyde	[$([CX3H][#6]),$([CX3H2])]=[OX1]
    Ketone	[#6][CX3](=[OX1])[#6]
    Thioaldehyde	[$([CX3H][#6]),$([CX3H2])]=[SX1]
    Thioketone	[#6][CX3](=[SX1])[#6]
    Imine	[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16])]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])]
    Carboxylic acid	[CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]
    Amide	[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]
    Nitro	[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]
    Nitrile	[NX1]#[CX2]
    Isonitrile	[CX1-]#[NX2+]
    Isocyanate	[NX2]=[CX2]=[OX1]
    Cyanate	[OX2][CX2]#[NX1]
    Isothiocyanate	[NX2]=[CX2]=[SX1]
    Thiocyanate	[SX2][CX2]#[NX1]
    Carbodiimide	[NX2]=[CX2]=[NX2]
    '''
    
    fl_gps = ['[CX4H3]',
              '[CX4H2]',
              '[CX4H]',
              '[#6]=[#6]',
              '[CX2]#[CX2]',
              '[CX3]=[CX2]=[CX3]',
              '[ClX1]',
              '[FX1]',
              '[BrX1]',
              '[IX1]',
              '[OX2H][c]',
              '[OX2H]',
              '[SX2H]',
              '[#6]-[#8]-[#6]',
              '[#6]-[#16]-[#6]',
              '[NX3;H2]',
              '[NX3;H]',
              '[SX2D2][SX2D2]',
              '[$([CX3H][#6]),$([CX3H2])]=[OX1]',
              '[#6][CX3](=[OX1])[#6]',
              '[$([CX3H][#6]),$([CX3H2])]=[SX1]',
              '[#6][CX3](=[SX1])[#6]',
              '[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16])]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])]',
              '[CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-])]',
              '[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]',
              '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',
              '[NX1]#[CX2]',
              '[CX1-]#[NX2+]',
              '[NX2]=[CX2]=[OX1]',
              '[OX2][CX2]#[NX1]',
              '[NX2]=[CX2]=[SX1]',
              '[SX2][CX2]#[NX1]',
              '[NX2]=[CX2]=[NX2]'
                ]
    return fl_gps




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
    

    
    
def search_substructures(input_path,output_path, cid_no, fl_gps):
    #This function serch for the presence of different functional groups/sub-structures for a given set of molecules in sdf format.
    file = output_path
    
    sub_strutures = fl_gps # ['[OX2H]', '[NX3;H2]', '[NX3;H]', '[SX2H]']
    header = ['cid_no']+sub_strutures
    cid = cid_no
    
    with open(file, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header of csv file
        writer.writerow(header)
        
        #open sdf file
        suppl = Chem.SDMolSupplier(input_path,removeHs = False)
        i=0
        for mol in suppl:   
        
            ful_gps = [cid[i]]
            for struc in sub_strutures:
                
                patt = Chem.MolFromSmarts(struc)
                presence = int(mol.HasSubstructMatch(patt))
                ful_gps.append(presence)
                
                
            print(ful_gps)
            writer.writerow(ful_gps)
            i = i+1
               
                    
                    
                    
                    
def main():
    input_path = input("Enter Path for input sdf file: ")
    output_path = input("Enter Path for output csv file: ")
    cid_no = retrieve_cid_identifiers(input_path)
    fl_gps = sub_structures()
    search_substructures(input_path,output_path, cid_no, fl_gps)
    
  


if __name__ == '__main__':
    main()
               
                    
                    




