# ----------------Molecular Weight Calculator-----------------------------------------------------------------
# Calculate molecular weight of a given molecular formula
# Required user inputs : 
#   No of elements present in the formula, Symbol of each elemnt and no of atoms from each element
#
# Author: Erandika Karunaratne
# Date : 2022-04-10
# ------------------------------------------------------------------------------------------------------------

import pandas as pd


#file =".....Atomic_weight.csv"
#df_mw = pd.read_csv(file)
#df_mw.dtypes
#df_mw.columns
#aw =pd.Series(df_mw.Atomic_Wt.values,index=df_mw.Symbol).to_dict()
#atomic_weight = aw

# ---Atomic Weights--------------

atomic_weight = {'Ac': 227.02775,
 'Ag': 107.8682,
 'Al': 26.981538399999994,
 'Am': 243.06137999999999,
 'Ar': 39.95,
 'As': 74.921595,
 'At': 209.98715,
 'Au': 196.96657,
 'B': 10.81,
 'Ba': 137.327,
 'Be': 9.0121831,
 'Bh': 270.133,
 'Bi': 208.9804,
 'Bk': 247.07031,
 'Br': 79.904,
 'C': 12.011,
 'Ca': 40.078,
 'Cd': 112.414,
 'Ce': 140.116,
 'Cf': 251.07959,
 'Cl': 35.45,
 'Cm': 247.07035,
 'Cn': 286.17900000000003,
 'Co': 58.93319399999999,
 'Cr': 51.9961,
 'Cs': 132.90545196,
 'Cu': 63.54600000000001,
 'Db': 268.126,
 'Ds': 282.166,
 'Dy': 162.5,
 'Er': 167.25900000000001,
 'Es': 252.083,
 'Eu': 151.964,
 'F': 18.998403,
 'Fe': 55.845,
 'Fl': 290.192,
 'Fm': 257.09511000000003,
 'Fr': 223.01972999999998,
 'Ga': 69.723,
 'Gd': 157.25,
 'Ge': 72.63,
 'H': 1.008,
 'He': 4.0026019999999995,
 'Hf': 178.486,
 'Hg': 200.592,
 'Ho': 164.930329,
 'Hs': 269.1336,
 'I': 126.90446999999999,
 'In': 114.818,
 'Ir': 192.217,
 'K': 39.0983,
 'Kr': 83.79799999999999,
 'La': 138.90547,
 'Li': 6.94,
 'Lr': 266.12,
 'Lu': 174.9668,
 'Lv': 293.205,
 'Mc': 290.19599999999997,
 'Md': 258.09843,
 'Mg': 24.305,
 'Mn': 54.93804300000001,
 'Mo': 95.95,
 'Mt': 277.154,
 'N': 14.007,
 'Na': 22.989769,
 'Nb': 92.90637,
 'Nd': 144.24200000000002,
 'Ne': 20.1797,
 'Nh': 286.182,
 'Ni': 58.6934,
 'No': 259.101,
 'Np': 237.04817200000002,
 'O': 15.999,
 'Og': 295.216,
 'Os': 190.23,
 'P': 30.973761,
 'Pa': 231.03588,
 'Pb': 207.2,
 'Pd': 106.42,
 'Pm': 144.91276000000002,
 'Po': 208.98243,
 'Pr': 140.90766000000002,
 'Pt': 195.084,
 'Pu': 244.0642,
 'Ra': 226.02541000000002,
 'Rb': 85.4678,
 'Re': 186.207,
 'Rf': 267.122,
 'Rg': 282.16900000000004,
 'Rh': 102.90549,
 'Rn': 222.01757999999998,
 'Ru': 101.07,
 'S': 32.06,
 'Sb': 121.76,
 'Sc': 44.955907,
 'Se': 78.971,
 'Sg': 269.128,
 'Si': 28.085,
 'Sm': 150.36,
 'Sn': 118.71,
 'Sr': 87.62,
 'Ta': 180.94788,
 'Tb': 158.925354,
 'Tc': 98.9062,
 'Te': 127.6,
 'Th': 232.0377,
 'Ti': 47.867,
 'Tl': 204.38,
 'Tm': 168.93421899999998,
 'Ts': 294.211,
 'U': 238.02891,
 'V': 50.9415,
 'W': 183.84,
 'Xe': 131.293,
 'Y': 88.905838,
 'Yb': 173.045,
 'Zn': 65.38,
 'Zr': 91.22399999999999}


def calculate_molecular_weight(no_atoms):
    """This function calculates and return molecular weight of Molecular formula"""
    molecular_formula =""
    molecular_formula_dict = {}
    #atomic_weight_total =[]
    molecular_weight = 0
    
    
    for i in range(1,no_atoms+1):
        element = input(f"Enter element {i}(H,C,N etc.):")
        no_atoms_of_element = int(input("Enter no of atoms: "))
        
        if no_atoms_of_element ==1:
            molecular_formula = molecular_formula + element
        else:
            molecular_formula = molecular_formula + element + str(no_atoms_of_element)
        
        molecular_formula_dict[element] = no_atoms_of_element
        individual_atomic_weight = no_atoms_of_element * atomic_weight[element]
        #atomic_weight_total.append(individual_atomic_weight)
        molecular_weight = molecular_weight + individual_atomic_weight
        i = i+1
        
    print(f"Molecular weight of {molecular_formula}: {molecular_weight} a.u")
    return molecular_weight
    return molecular_formula



def main():
    no_atoms = int(input("Enter no of elements present in your molecular formula: "))
    calculate_molecular_weight(no_atoms)
    
    
if __name__ == "__main__":
    main()



