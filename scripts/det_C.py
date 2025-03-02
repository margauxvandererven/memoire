import matplotlib.pyplot as plt
from tueplots import fonts

from readmultispec import * #permet de lire les fichiers fits de spectre à échelle
from wavelen_work import * #permet de traiter les spectres synthétiques et de zoomer sur différentes parties du spectres
from zoom_raies import * #permet de zoomer sur les raies demandées avec les spectres synthétiques demandés
from minimisation_chi_2 import *
from stardata_BD22 import *
plt.rcParams.update(fonts.neurips2021())

element="CO"
raie_element={}

raie_propre = {
        15581.85 :[None, None],

}

nom="CO_/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14700-18000_BD-221742_"
ABU = [7.85, 7.9, 7.95, 8.0, 8.1]

synth = {   # "toutsans_OH":"sans",
            nom+"Cabu_7.85"+".conv":"log $\epsilon_C$ = 7.85",
            nom+"Cabu_7.9"+".conv":"log $\epsilon_C$ = 7.9",
            nom+"Cabu_7.95"+".conv":"log $\epsilon_C$ = 7.95",
            nom+"Cabu_8.0"+".conv":"log $\epsilon_C$ = 8.0",
            nom+"Cabu_8.1"+".conv":"log $\epsilon_C$ = 8.1",
                }

# synth_plot = { "toutsans_OH":"sans OH",
#             nom+"Oabu_8.56"+".conv":"log $\epsilon_O$ = 8.56",
#             # nom+"Oabu_8.58"+".conv":"log $\epsilon_O$ = 8.58",
#             nom+"Oabu_8.6"+".conv":"log $\epsilon_O$ = 8.6",
#             # nom+"Oabu_8.62"+".conv":"log $\epsilon_O$ = 8.62",
#             nom+"Oabu_8.65"+".conv":"log $\epsilon_O$ = 8.65",
#                 }

for wavelength in raie_propre:
    chi_plot_ABU(
            ABU, 
            chi_2(path, synth, stardata, wavelength, lines_BD22,
                name=element,start=raie_propre.get(wavelength)[0],end=raie_propre.get(wavelength)[1], plot=True)["chi_squared_values"], "C",wavelength, element, raie_element)
#     chi_2(path, synth_plot, stardata, wavelength, lines_BD22,
#                 name=element,start=raie_propre.get(wavelength)[0],end=raie_propre.get(wavelength)[1], plot=True)

# abu_plot(raie_element,"O")



# fausses_raies={"x":list(np.arange(14700,15500,20))}
# zoom_lines(fausses_raies, path, {"CO_/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_2_Cabu_8.44_ratio_90.conv":"C=8.44, ratio 90",
#                                  "CO_/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_2_Cabu_8.44_ratio_4.conv":"C=8.44, ratio 4"
#                                  },stardata,10,lines_BD22)


# raies={"13C17O" : [14742,
#                 #    14864.1, 14879.5, 14884.2,14921.2,14994,15208,15229,15304,
# 15442.3,15663]
# }
# zoom_lines(raies, path, {"CO_/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_2_Cabu_8.44_ratio_90.conv":"C=8.44, ratio 90",
#                         "CO_/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_2_Cabu_8.44_ratio_4.conv":"C=8.44, ratio 4",
#                         "CO_/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14700-15700_BD-221742_Cabu_8.44_ratio_70.conv":"C=8.44, ratio 70"
# }, stardata,5,lines_BD22)

