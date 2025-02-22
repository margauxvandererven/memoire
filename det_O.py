import matplotlib.pyplot as plt
from tueplots import fonts

from readmultispec import * #permet de lire les fichiers fits de spectre à échelle
from wavelen_work import * #permet de traiter les spectres synthétiques et de zoomer sur différentes parties du spectres
from zoom_raies import * #permet de zoomer sur les raies demandées avec les spectres synthétiques demandés
from minimisation_chi_2 import *
from stardata_BD22 import *
plt.rcParams.update(fonts.neurips2021())

element="OH"
raie_element={}

raie_propre = {
        15002.17: [15001.75,15002.7], 
        15391.12 : [15390.56, 15391.7], 
        15651.897 : [15651.35,15652.45], 
        # 16247.884 : [16247.65,16248.4],
        # 15236.75 : [15236.3,15237.3],
        # 16581.269 : [16580.74,16581.64] ,
        15003.15: [15002.67, 15003.7],
        15266.2: [15265.7,15266.65],
        17772.7: [17772.1, 17773.25],
        17423.85: [17423.39, 17424.35], 
        16662.21: [16661.9,16662.76], 
        16448.1: [16447.5, 16448.59],
        16368.15: [16367.6, 16368.7]
}

nom="OH/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_2_"
ABU = [8.56,8.58, 8.6,8.62,8.65]

synth = {   # "toutsans_OH":"sans",
            nom+"Oabu_8.56"+".conv":"log $\epsilon_O$ = 8.56",
            nom+"Oabu_8.58"+".conv":"log $\epsilon_O$ = 8.58",
            nom+"Oabu_8.6"+".conv":"log $\epsilon_O$ = 8.6",
            nom+"Oabu_8.62"+".conv":"log $\epsilon_O$ = 8.62",
            nom+"Oabu_8.65"+".conv":"log $\epsilon_O$ = 8.65",
                }

synth_plot = { "toutsans_OH":"sans OH",
            nom+"Oabu_8.56"+".conv":"log $\epsilon_O$ = 8.56",
            # nom+"Oabu_8.58"+".conv":"log $\epsilon_O$ = 8.58",
            nom+"Oabu_8.6"+".conv":"log $\epsilon_O$ = 8.6",
            # nom+"Oabu_8.62"+".conv":"log $\epsilon_O$ = 8.62",
            nom+"Oabu_8.65"+".conv":"log $\epsilon_O$ = 8.65",
                }

for wavelength in raie_propre:
    chi_plot_ABU(
            ABU, 
            chi_2(path, synth, stardata, wavelength, lines_BD22,
                name=element,start=raie_propre.get(wavelength)[0],end=raie_propre.get(wavelength)[1])["chi_squared_values"], "O",wavelength, element, raie_element)
    chi_2(path, synth_plot, stardata, wavelength, lines_BD22,
                name=element,start=raie_propre.get(wavelength)[0],end=raie_propre.get(wavelength)[1], plot=True)

abu_plot(raie_element,"O")