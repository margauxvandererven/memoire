from imports import *
import pync
# from utilities import *
from tueplots import fonts
from scipy import stats
plt.rcParams.update(fonts.neurips2021())


# filename = "Fe_lines_"+stardata.get("starname")+".txt"
filename = "Fe_final.txt"
filename2 = "Fe_lines_"+stardata.get("starname")+"_bacchus_4307.txt"
filename3 = "Fe_lines_"+stardata.get("starname")+"_bacchus_4258.txt"

with open("../results/"+filename, "r") as fichier:
    chi_final_data_4000 = json.load(fichier)
with open("../results/"+filename2, "r") as fichier:
    chi_final_data2 = json.load(fichier)
with open("../results/"+filename3, "r") as fichier:
    chi_final_data3 = json.load(fichier)

with open("Fe_final_4125.txt", "r") as fichier:
    chi_final_data_4125 = json.load(fichier)

with open("Fe_final_3875.txt", "r") as fichier:
    chi_final_data_3875 = json.load(fichier)

with open("Fe_final_3800.txt", "r") as fichier:
    chi_final_data_3800 = json.load(fichier)

with open("../data_lines/raies.txt", "r") as fichier:
    config=json.load(fichier)
raie_Fe=config["Fe_lines"]
# print(raie_Fe)

# print(raie_Fe)

# abu_plot(chi_final_data,"Fe",size_police=20,
#      save="Fe_abu_4000"
#     )

# abu_plot(chi_final_data2,"Fe",size_police=20,
#      save="Fe_abu_bacchus_4307"
#     )
# abu_plot(chi_final_data3,"Fe",size_police=20,
#      save="Fe_abu_bacchus_4258"
#     )

def plot_Teff(lines_data, chi_final_data, size_police=12, save=None, chi_exc_range=None):
    """
    lines_data : {line:[start, end, exc pot, log gf]}
    chi_final_data : {line:[start, end, abundance,error_abu_minus,error_abu_plus,chi2_red,p_value]}
    """
    element_abu="Fe"
    x_vals = []
    y_vals = []
    errors_minus = [] 
    errors_plus = []

    for i in chi_final_data:
        val=chi_final_data.get(i)
        p_value=np.float64(val[-1])
        if 0.05<p_value<0.95 and val[3]>0 and val[4]>0:
            exc_pot = lines_data.get(i)[2]
            if chi_exc_range is None or (chi_exc_range[0] <= exc_pot <= chi_exc_range[1]):
                x_vals.append(exc_pot)
                y_vals.append(val[2])
                errors_minus.append(0.07)  
                errors_plus.append(0.07) 

    x_vals= np.array(x_vals)
    errors = np.array([errors_minus, errors_plus])
    y_vals=np.array(y_vals)

    f = plt.figure(figsize=(12, 5))
    gs = f.add_gridspec(1)
    ax = gs.subplots(sharex=False, sharey=True)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
    confidence = 0.95
    degrees_of_freedom = len(x_vals) - 2
    t_value = stats.t.ppf((1 + confidence) / 2, degrees_of_freedom)
    # Calcul des intervalles de confiance pour la pente et l'ordonnée
    slope_ci = t_value * std_err

    # Calcul de la ligne de régression et des intervalles de confiance
    x_line = np.linspace(min(x_vals), max(x_vals), 100)
    y_line = slope * x_line + intercept
    
    # Calcul de l'erreur standard de la régression
    y_pred = slope * x_vals + intercept
    residuals = y_vals - y_pred
    std_residuals = np.std(residuals)
    
    # Test de normalité des résidus
    _, normality_p_value = stats.normaltest(residuals)
    
    # Calcul du R² ajusté
    r_squared = r_value**2
    n = len(x_vals)
    p = 1  # nombre de variables explicatives
    r_squared_adj = 1 - (1 - r_squared) * (n - 1) / (n - p - 1)

    ax.errorbar(x_vals, y_vals,yerr=errors, color="darkblue", label=f"IR : 4000 K", capsize=5,  fmt='o')

    ax.plot(x_line, y_line, color='firebrick', 
            label=f'Régression: y = ({slope:.3f}±{slope_ci:.3f})x + {intercept:.3f}')
    
    # Ajout des statistiques sur le graphique
    stats_text = f'R² ajusté = {r_squared_adj:.3f}\n'
    stats_text += f'p-value = {p_value:.3e}\n'
    stats_text += f'Erreur std = {std_residuals:.3f}\n'
    stats_text += f'Normalité p-value = {normality_p_value:.3e}'
    ax.hlines(np.mean(y_vals), min(x_vals), max(x_vals), color='chocolate', linestyle='--', label="Abondance moyenne")
    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel("$\chi_{exc}$ [eV]", fontsize=size_police)
    ax.set_ylabel(f"$\\log{{\\epsilon_{{{element_abu}}}}}$", fontsize=size_police)

    ax.xaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    ax.xaxis.set_tick_params(direction='in', length=5, which='minor', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    confidence_band = stats.t.ppf(0.975, len(x_vals)-2) * std_residuals
    ax.fill_between(x_line, y_line - confidence_band, y_line + confidence_band, 
                color='gray', alpha=0.1)
    ax.tick_params(axis = 'both', labelsize = size_police)

    if save:
        plt.savefig(save+".pdf", dpi=600, bbox_inches='tight', transparent=True)
    plt.show()

    # Afficher les résultats détaillés
    print("\nRésultats de la régression linéaire:")
    print(f"Pente = {slope:.3f} ± {slope_ci:.3f}")
    print(f"Ordonnée à l'origine = {intercept:.3f}")
    print(f"R² ajusté = {r_squared_adj:.3f}")
    print(f"P-value = {p_value:.3e}")
    print(f"Erreur standard des résidus = {std_residuals:.3f}")
    print(f"Test de normalité des résidus p-value = {normality_p_value:.3e}")
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared_adj': r_squared_adj,
        'p_value': p_value,
        'std_residuals': std_residuals,
        'normality_p_value': normality_p_value
    }

plot_Teff(raie_Fe, chi_final_data_3875, size_police=18, save="../output/Teff_3875")

# plot_Teff(raie_Fe, chi_final_data, size_police=18, 
#           save="../présentation/images/comparaison_modèles", chi_exc_range=(4.5, 7)
#           )

# data = get_ew_atom(ew_limit=1e-10, Teff=4000, particular_element="Fe I")["data"]
# excitation_dict = {item[0]: [item[1], item[3]] for item in data}
# raies_Fe = {wl: excitation_dict.get(wl, "Non trouvé") for wl in raie_tot_Fe}


# model="4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod"
# model="s4125_g+1.0_m2.0_t02_st_z-0.30_a+0.12_c+0.00_n+0.00_o+0.12_r+0.00_s+0.00.int"
# model='s3875_g+1.0_m2.0_t02_st_z-0.30_a+0.12_c+0.00_n+0.00_o+0.12_r+0.00_s+0.00.int'
model='s3800_g+1.0_m2.0_t02_st_z-0.30_a+0.12_c+0.00_n+0.00_o+0.12_r+0.00_s+0.00.int'
# model="4258g2.04m1.0z-0.45_BD221742.int"
# model="4307g2.29m1.0z-0.37_BD221742.int"

ABU=[6.9, 6.85, 7.35, 7.3, 7.4, 7.5, 7.0, 7.2, 7.1, 7.25, 7.7, 7.24, 7.23, 7.22, 7.21, 7.19]
# analyse_chi2(raie_Fe, ABU, "Fe", "final", stardata, lines_BD22,
#              minimisation=True,
#              abu_to_plot=[7.2],
#              name="Fe", 
#              save="../results/Fe_final", plot=True
#              )
# abu_to_plot=[7.2]

def analyse_abu_Fe(raies,save=None, abu_to_plot=None, minimisation=None, plot=None):

    chi_final={}
    variable="Fe"
    spectral_lines=lines_BD22
    name="log $\epsilon_{Fe}$"
    path_to_synth2="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
    path_to_synth="/Users/margauxvandererven/Unif/memoire_local/Turbospectrum_NLTE-master/COM/syntspec/synth_margaux/"
    synth={}
    range="14630-22900"

    for abu in ABU:
        synth[model+"_"+range+"_"+variable+"abu_"+str(abu)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
    # print(synth)

    for wavelength in raies:
        print(f"Processing line {wavelength}")
        wavelength=np.float64(wavelength)
        
        if minimisation is not None:
            chi_squared_values = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])["chi_squared_values"]
            dof = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])["dof"]
            chi_minimisation_ABU(ABU, chi_squared_values, variable,wavelength, name, chi_final,dof=dof, plot=plot, 
                                     save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final_Fe_3800/"+str(wavelength)+"/"+str(wavelength)+"_minimisation"
                                )
        if abu_to_plot is not None:
                synth_plot={}
                for abu2 in abu_to_plot:
                    synth_plot[model+"_"+range+"_"+variable+"abu_"+str(abu2)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                    if wavelength < 18500:
                        synth_plot["../../../../../../../../Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/BD-221742b/Fe/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_Feabu_-20.conv"]= "sans Fe"
                    else:
                        synth_plot["../../../../../../../../Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/BD-221742b/Fe/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_19750-22900_BD-221742_Feabu_-20.conv"]= "sans Fe"
                
                start=raies.get(str(wavelength))[0]
                end=raies.get(str(wavelength))[1]       
                plot_zone_chi2_simple(wavelength, path_to_synth, synth_plot, stardata, spectral_lines, size_police=14,size_trace=(1., 8),name=name, start=start,end=end,
                                    # save="/Users/margauxvandererven/OneDrive  - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/"+str(wavelength)+"/"+str(wavelength)+"_zone", plot=plot
                                    save="/Users/margauxvandererven/OneDrive  - Université Libre de Bruxelles/memoire/output/final_Fe_3800/"+str(wavelength)+"/"+str(wavelength)+"_zone", plot=plot
                                    )
    keys_to_remove = []
    for k, v in chi_final.items():
        if len(v) != 8:
            keys_to_remove.append(k)
    
    for k in keys_to_remove:
        del chi_final[k]

    if minimisation is not None:
        abu_plot(chi_final,variable,size_police=20,
                #  save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/Oabu_"+round,
                    save="../rédaction/images/plot_abu/Fe_final_3800.pdf")

    if save:
        with open(save+".txt", "w") as fichier:
            json.dump(chi_final, fichier, indent=4, ensure_ascii=False)

        df = pd.DataFrame({
        'wavelength_center': list(chi_final.keys()),
        'wavelength_start': [val[0] for val in chi_final.values()],
        'wavelength_end': [val[1] for val in chi_final.values()],
        'abundance': [val[2] for val in chi_final.values()],
        'error_abu_minus': [val[3] for val in chi_final.values()],
        'error_abu_plus': [val[4] for val in chi_final.values()],
        'chi2_red': [val[5] for val in chi_final.values()],
        # '$R^2$': [val[6] for val in chi_final.values()],
        'p_value': [val[7] for val in chi_final.values()],
        })
        df.to_csv(save+".csv", index=False)

# print(raie_Fe.get(str(22493.67)))

# analyse_abu_Fe(raie_Fe,save="Fe_final_3800", abu_to_plot=[7.2],minimisation=True)


# print(raie_Fe.get(str(17012.729)))
# analyse_abu({"22493.67": raie_Fe.get(str(22493.67))}, abu_to_plot=abu_to_plot)

# abu_plot(chi_final_data_3800,"Fe",size_police=20,
#                 #  save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/Oabu_"+round,
#                     save="../rédaction/images/plot_abu/Fe_final_3800.pdf")

# Afficher une notification avec pync
pync.notify("Votre script a terminé son exécution.", title="Script terminé")