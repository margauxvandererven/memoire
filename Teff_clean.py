from imports import *

filename = "Fe_lines_"+stardata.get("starname")+"_bacchus_4258.txt"

# with open(filename, "r") as fichier:
#     chi_final_data = json.load(fichier)

with open("raies.txt", "r") as fichier:
    config=json.load(fichier)

raie_Fe=config["Fe_lines"]
print(raie_Fe)
# abu_plot(chi_final_data,"Fe",size_police=20,
#     #  save=variable+"_abu"
#     )
# plot_Teff(lines_data, chi_final_data, size_police=12, save=None)

# data = get_ew_atom(ew_limit=1e-10, Teff=4000, particular_element="Fe I")["data"]
# excitation_dict = {item[0]: [item[1], item[3]] for item in data}
# raies_Fe = {wl: excitation_dict.get(wl, "Non trouvé") for wl in raie_tot_Fe}


# model="4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod"
model="4258g2.04m1.0z-0.45_BD221742.int"
ABU=[7.35, 7.3, 7.4, 7.5, 7.0, 7.2, 6.9, 7.1, 7.25, 7.7]
abu_to_plot=[7.2, 7.3]

def analyse_abu(raies,save=None, abu_to_plot=None):
    chi_final={}
    variable="Fe"
    spectral_lines=lines_BD22
    name="log $\epsilon_{Fe}$"
    path_to_synth="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"

    for wavelength in raies:
        wavelength=np.float64(wavelength)
        range="14630-22900"
        synth={}
        for abu in ABU:
            synth[model+"_"+range+"_"+variable+"abu_"+str(abu)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
        if abu_to_plot is not None:
            synth_plot={}
            for abu2 in abu_to_plot:
                synth_plot["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_"+range+"_"+variable+"abu_"+str(abu2)+"_"+round+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                synth_plot["../../syntspec/BD-221742b/sans_CO_mol2"]= "sans "+name
            plot_zone_chi2(wavelength, path_to_synth, synth_plot, stardata, spectral_lines,axes=(True,True), size_police=20,size_trace=(1.8, 10),name=name, start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])
        chi_squared_values = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])["chi_squared_values"]
        chi_minimisation_ABU(ABU, chi_squared_values, variable,wavelength, name, chi_final, plot=True)
    if save:
        with open(save, "w") as fichier:
            json.dump(chi_final, fichier, indent=4, ensure_ascii=False)

# analyse_abu(raie_Fe, filename)