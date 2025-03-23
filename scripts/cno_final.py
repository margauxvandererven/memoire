from imports import *
import pync
repertory_memoire="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/"
# path_to_synth="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/BD-221742b/"

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
    # 8.54, 
    8.55, 
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
            #   save="../results/OH_final_"+round+".txt"
              )


# raies_OH_final={14661.144: [14660.6, 14661.59, -5.989],
#  15002.153: [15001.67, 15002.67, -5.653],
#  15003.12: [15002.67, 15003.64, -5.653],
#  15130.921: [15130.32, 15131.53, -5.572],
#  15266.168: [15265.54, 15266.76, -5.5],
#  15278.525: [15278.03, 15279.16, -5.453],
#  15391.205: [15390.49, 15391.93, -5.512],
#  15409.17: [15408.47, 15409.77, -5.435],
#  15428.401: [15427.96, 15429.06, -5.425],
#  15429.688: [15429.03, 15430.13, -5.153],
#  15505.746: [15504.8, 15506.66, -5.378],
#  15568.782: [15568.17, 15569.43, -5.337],
#  15651.897: [15651.33, 15652.41, -5.203],
#  15719.695: [15719.09, 15720.36, -5.32],
#  15755.52: [15754.86, 15756.04, -5.175],
#  15756.53: [15756.06, 15757.1, -5.175],
#  16052.766: [16052.15, 16053.29, -4.976],
#  16247.884: [16247.5, 16248.78, -5.177],
#  16312.92: [16311.8, 16313.42, -5.077],
#  16347.493: [16347.12, 16348.07, -5.002],
#  16368.136: [16367.51, 16368.76, -4.858],
#  16662.2: [16661.77, 16662.85, -5.069],
#  16729.784: [16729.47, 16730.32, -4.793],
#  16904.278: [16903.62, 16904.82, -4.712],
#  17322.25: [17321.83, 17322.84, -4.628],
#  17423.859: [17423.35, 17424.43, -4.503]}

# nom="OH/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_2_"
# ABU = [8.56, 8.58, 8.6, 8.62, 8.65]
# synth_OH_1={}
# for abu in ABU : 
#     synth_OH_1[nom+"Oabu_"+str(abu)+".conv"]="log $\epsilon_O$ = "+str(abu)
# chi_final_OH_1={}
# for wavelength in raies_OH_final:
#     start=raies_OH_final.get(wavelength)[0]
#     end=raies_OH_final.get(wavelength)[1]
#     chi_squared_values = chi_2(path_to_synth, synth_OH_1, stardata, wavelength, lines_BD22,chi_final_OH_1,name="16OH",start=start,end=end)["chi_squared_values"]
#     dof=chi_2(path_to_synth, synth_OH_1, stardata, wavelength, lines_BD22,chi_final_OH_1,name="16OH",start=start,end=end)["dof"]
#     chi_minimisation_ABU(ABU, chi_squared_values, "O",wavelength, "16OH", chi_final_OH_1,dof=dof, plot=True, 
#                         #  save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/"+str(wavelength)+"/"+str(wavelength)+"_minimisation"
#                          )    
#     # synth_plot_OH_1=np.copy(synth_OH_1)
#     # synth_plot_OH_1["toutsans_OH"]="sans OH"
#     # plot_zone_chi2_simple(wavelength, path_to_synth, synth_plot_OH_1, stardata, spectral_lines=lines_BD22, size_police=14,size_trace=(1., 8),name="16OH", start=start,end=end,
#     #                                 # save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/"+str(wavelength)+"/"+str(wavelength)+"_zone"
#     #                                 )
# abu_plot(chi_final_OH_1,"C",size_police=20)