from imports import *
from minimisation_chi_2 import *

def calculate_C_abu(iteration, ABU):
    with open('../data_lines/raies_moleculaires_test.txt', 'r') as f:
        data = json.load(f)
    raies=data["12C16O_lines"]

    analyse_chi2(raies, ABU, "C",stardata,lines_BD22, minimisation=True, abu_to_plot=None,plot=None, name=None,save=f"chi_CO_{iteration}")

    with open(f"chi_CO_{iteration}.txt", "r") as fichier:
        raies_element = json.load(fichier)

    y_vals=[]
    for key in raies_element.keys():
        val=raies_element[key]
        p_value=np.float64(val[-1])
        if 0.05<p_value<0.95 and val[3]>0 and val[4]>0:
            y_vals.append(val[2]) 

    return np.mean(y_vals) 