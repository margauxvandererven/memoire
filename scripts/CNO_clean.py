from imports import *
import pync
repertory_memoire="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/"

#Abu de litt.
# C = 0.35+8.39-0.3 #asp2007
# N=-0.1+7.78-0.3 #asp2007
# print(C,N)

with open("../data_lines/raies_moleculaires.txt", "r") as fichier:
            config= json.load(fichier)

# -------------------------------------------------------------------------------------
# Oxygène
# -------------------------------------------------------------------------------------
OH_raies=config["OH_lines"]

O_ABU = [
    8.54, 8.55, 
    8.56, 
    8.57, 
    8.58,
    8.59, 
    8.60,
    8.61,
    8.62,
    8.63,
    8.64,
    8.65
    ]
round="first"

analyse_chi2(OH_raies, O_ABU, "O", round,stardata,lines_BD22,
             minimisation=True,
             abu_to_plot=[8.57, 8.60,8.63],
             name="OH", 
              save="../results/OH_final_"+round+".txt"
              )

# with open("../results/OH_"+round+".txt", "r") as fichier:
#         chi_final_OH=json.load(fichier)

# -------------------------------------------------------------------------------------
# Carbone
# -------------------------------------------------------------------------------------
raies_12CO = config["12C16O_lines"]

C_ABU=[7.5, 7.85, 7.9, 7.95, 8.0, 8.1, 8.2, 8.3, 8.4, 7.7, 7.8, 7.65]

# analyse_chi2_CO(raies_12CO,C_ABU, "C", "first", stardata, lines_BD22,minimisation=True,abu_to_plot=[7.9, 8.0, 7.8],name="12CO",save="../results/CO_first.txt")



# pync.notify("Votre script a terminé son exécution.", title="Script terminé")