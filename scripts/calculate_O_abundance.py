from imports import *
from minimisation_chi_2 import *
import sys


def calculate_O_abu(iteration, abundances_str):
    ABU = [float(x) for x in abundances_str.split(',')]
    with open('../data_lines/raies_moleculaires_test.txt', 'r') as f:
        data = json.load(f)
    raies=data["OH_lines"]

    analyse_chi2(raies, ABU, "O",stardata,lines_BD22,minimisation=True, abu_to_plot=None,plot=None, name="OH",save=f"chi_OH_{iteration}")

    with open(f"chi_OH_{iteration}.txt", "r") as fichier:
        raies_element = json.load(fichier)

    y_vals=[]
    for key in raies_element.keys():
        val=raies_element[key]
        p_value=np.float64(val[-1])
        if 0.05<p_value<0.95 and val[3]>0 and val[4]>0:
            y_vals.append(val[2]) 

    return np.mean(y_vals) 

if __name__ == "__main__":
    # Vérification du nombre d'arguments
    if len(sys.argv) != 3:
        print("Usage: python calculate_O_abundance.py <iteration> <abundances>")
        print("Example: python calculate_O_abundance.py 1 '8.1,8.2,8.3,8.4,8.5'")
        sys.exit(1)
    
    # Récupération des arguments
    iteration = sys.argv[1]
    abundances = sys.argv[2]
    
    # Calcul et affichage du résultat
    result = calculate_O_abu(iteration, abundances)
    print(result)

