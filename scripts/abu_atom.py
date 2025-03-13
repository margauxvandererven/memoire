from imports import *
import pync
from tueplots import fonts
plt.rcParams.update(fonts.neurips2021())


with open("../data_lines/raies.txt", "r") as fichier:
    config=json.load(fichier)

raie_Ti=config["Ti1_lines"]
raie_Ti_ion=config["Ti2_lines"]
print(raie_Ti)
model="4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod"


ABU=[5, 4.9, 4.7, 4.5, 4, 4.2]
def analyse_abu_Ti(raies,save=None, abu_to_plot=None, minimisation=None):
    chi_final={}
    variable="Ti"
    spectral_lines=lines_BD22
    name="log $\epsilon_{Ti}$"
    path_to_synth="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"

    for wavelength in raies:
        print(wavelength)
        wavelength=np.float64(wavelength)
        # range="14820-23200"
        range=str(wavelength-5)+"-"+str(wavelength+5)
        if minimisation is not None:
            synth={}
            for abu in ABU:
                synth[model+"_"+range+"_"+variable+"abu_"+str(abu)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
            chi_squared_values = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])["chi_squared_values"]
            chi_minimisation_ABU(ABU, chi_squared_values, variable,wavelength, name, chi_final,
                                  plot=True
                                )
        if abu_to_plot is not None:
            synth_plot={}
            for abu2 in abu_to_plot:
                synth_plot[model+"_"+range+"_"+variable+"abu_"+str(abu2)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                synth_plot["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14820-23200_Tiabu_-20.conv"]= "sans Ti"
            plot_zone_chi2(wavelength, path_to_synth, synth_plot, stardata, spectral_lines,axes=(True,True), size_police=20,size_trace=(1.8, 10),name=name, start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])
            
        if save:
            with open(save, "w") as fichier:
                json.dump(chi_final, fichier, indent=4, ensure_ascii=False)

# raie_particulière="22890.01" #manque molécule ?TiI
raie_particulière="14711.708"

analyse_abu_Ti({raie_particulière:raie_Ti_ion.get(raie_particulière)}, abu_to_plot=[4.7, 4.5, 4.2],
            #    minimisation=True, save="../results/Ti.txt"
               )

# "Ti2_lines": {
#         "14631.66": [
#             14631.086,
#             14632.325,
#             3.124,
#             -1.724
#         ],
#         "14711.708": [
#             14711.108,
#             14712.379,
#             3.095,
#             -1.855
#         ],
#         "15231.203": [
#             15230.79,
#             15231.48,
#             3.123,
#             -3.617
#         ],
#         "15873.838": [
#             15873.141,
#             15874.606,
#             3.123,
#             -2.061
#         ],
#         "16621.818": [
#             16621.399,
#             16622.27,
#             3.123,
#             -2.567
#         ],
#         "17015.951": [
#             17015.751,
#             17016.241,
#             3.095,
#             -3.353
#         ]

analyse_abu_Ti(raie_Ti_ion, abu_to_plot=[4.7, 4.5, 4.2],
            #    minimisation=True, save="../results/Ti_ion.txt"
               )

# with open("../results/Ti_ion.txt", "r") as fichier:
#     chi_final_data = json.load(fichier)

# abu_plot(chi_final_data,"Ti",size_police=20,
#     #  save="Fe_abu_4000"
#     )