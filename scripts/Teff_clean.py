from imports import *
import pync
# from utilities import *
from tueplots import fonts
from scipy import stats
plt.rcParams.update(fonts.neurips2021())


filename = "Fe_lines_"+stardata.get("starname")+".txt"
filename2 = "Fe_lines_"+stardata.get("starname")+"_bacchus_4307.txt"
filename3 = "Fe_lines_"+stardata.get("starname")+"_bacchus_4258.txt"

with open("../results/"+filename, "r") as fichier:
    chi_final_data = json.load(fichier)
with open("../results/"+filename2, "r") as fichier:
    chi_final_data2 = json.load(fichier)
with open("../results/"+filename3, "r") as fichier:
    chi_final_data3 = json.load(fichier)


with open("../data_lines/raies.txt", "r") as fichier:
    config=json.load(fichier)

raie_Fe=config["Fe_lines"]
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

def plot_Teff(lines_data, chi_final_data, size_police=12, save=None):
    """
    lines_data : {line:[start, end, exc pot, log gf]}
    chi_final_data : {line:[start, end, log $\epsilon$, $\chi2$]}
    """
    element_abu="Fe"
    x_vals = []
    y_vals = []
    y_vals2 = []
    y_vals3 = []
    for i in lines_data:
        x_vals.append(lines_data.get(i)[2])
        y_vals.append(chi_final_data.get(i)[2])
        y_vals2.append(chi_final_data2.get(i)[2])
        y_vals3.append(chi_final_data3.get(i)[2])

    x_vals= np.array(x_vals)
    errors=np.ones(len(y_vals))*0.1
    y_vals=np.array(y_vals)
    y_vals2=np.array(y_vals2)
    y_vals3=np.array(y_vals3)

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

    ax.errorbar(x_vals, y_vals,yerr=errors, color="gray", label=f"IR : 4000 K", capsize=5,  fmt='o')
    # ax.scatter(x_vals, y_vals2, color="gray", label=f"Vis : 4307 K",marker='x', linewidths=0.5)
    # ax.scatter(x_vals, y_vals3, color="gray", label=f"Vis : 4258 K",marker='+', linewidths=1.2)

    ax.plot(x_line, y_line, color='darkblue', 
            label=f'Régression: y = ({slope:.3f}±{slope_ci:.3f})x + {intercept:.3f}')
    
    # Ajout des statistiques sur le graphique
    stats_text = f'R² ajusté = {r_squared_adj:.3f}\n'
    stats_text += f'p-value = {p_value:.3e}\n'
    stats_text += f'Erreur std = {std_residuals:.3f}\n'
    stats_text += f'Normalité p-value = {normality_p_value:.3e}'
    ax.hlines(7.2, min(x_vals), max(x_vals), color='indianred', linestyle='--', label="Abondance moyenne")
    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel("$\chi_{exc}$ [eV]", fontsize=size_police)
    ax.set_ylabel(f"$\\log{{\\epsilon_{{{element_abu}}}}}$", fontsize=size_police)

    ax.xaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    ax.xaxis.set_tick_params(direction='in', length=5, which='minor', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)

    # model = LinearRegression()
    # model.fit(x_vals, y_vals)
    # X_line = np.linspace(min(x_vals), max(x_vals), 100).reshape(-1, 1) 
    # Y_line = model.predict(X_line) 
    # slope = model.coef_[0]  
    # intercept = model.intercept_

    # coeffs, cov = np.polyfit(x_vals, y_vals, 1, cov=True)
    # m, b = coeffs  # Pente et ordonnée à l'origine
    # err_m, err_b = np.sqrt(np.diag(cov)) 
    # X_line = np.linspace(min(x_vals), max(x_vals), 100).reshape(-1, 1) 
    # Y_line = m * X_line + b 

    # ax.plot(X_line, Y_line, color='darkblue', label=f"Régression: y = {m:.2f}x + {b:.2f}")

    # plt.fill_between(X_line, (m - err_m) * X_line + (b - err_b), 
    #                       (m + err_m) * X_line + (b + err_b), 
    #                       color='red', alpha=0.2, label="Incertitude sur la régression")

    # model2 = LinearRegression()
    # model2.fit(x_vals, y_vals2)
    # Y_line2 = model2.predict(X_line) 
    # slope2 = model2.coef_[0]  
    # intercept2 = model2.intercept_ 
    # ax.plot(X_line, Y_line2, color='brown', label=f'Vis : 4307 K - 2.29')

    # model3 = LinearRegression()
    # model3.fit(x_vals, y_vals3)
    # Y_line3 = model3.predict(X_line) 
    # slope3 = model3.coef_[0]  
    # intercept3 = model3.intercept_ 
    # ax.plot(X_line, Y_line3, color='orange', label=f'Vis : 4258 K - 2.04')

    ax.tick_params(axis = 'both', labelsize = size_police)

    # plt.legend(ncol=2, framealpha=0.2, fontsize=size_police)
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
    # return f'y = {slope:.2f}x + {intercept:.2f}'


plot_Teff(raie_Fe, chi_final_data, size_police=18, 
          save="../présentation/images/comparaison_modèles"
          )

# data = get_ew_atom(ew_limit=1e-10, Teff=4000, particular_element="Fe I")["data"]
# excitation_dict = {item[0]: [item[1], item[3]] for item in data}
# raies_Fe = {wl: excitation_dict.get(wl, "Non trouvé") for wl in raie_tot_Fe}


model="4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod"
# model="4258g2.04m1.0z-0.45_BD221742.int"
# model="4307g2.29m1.0z-0.37_BD221742.int"

# ABU=[7.35, 7.3, 7.4, 7.5, 7.0, 7.2, 6.9, 7.1, 7.25, 7.7, 7.6]
# abu_to_plot=[7.2]

def analyse_abu_Fe(raies,save=None, abu_to_plot=None, minimisation=None):
    chi_final={}
    variable="Fe"
    spectral_lines=lines_BD22
    name="log $\epsilon_{Fe}$"
    path_to_synth="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"

    for wavelength in raies:
        wavelength=np.float64(wavelength)
        range="14630-22900"
        if minimisation is not None:
            synth={}
            for abu in ABU:
                synth[model+"_"+range+"_"+variable+"abu_"+str(abu)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
            chi_squared_values = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])["chi_squared_values"]
            chi_minimisation_ABU(ABU, chi_squared_values, variable,wavelength, name, chi_final,
                                #   plot=True
                                )
    if abu_to_plot is not None:
            synth_plot={}
            for abu2 in abu_to_plot:
                synth_plot[model+"_"+range+"_"+variable+"abu_"+str(abu2)+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                if wavelength < 18500:
                    synth_plot["../../syntspec/BD-221742b/Fe/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_Feabu_-20.conv"]= "sans Fe"
                else:
                    synth_plot["../../syntspec/BD-221742b/Fe/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_19750-22900_BD-221742_Feabu_-20.conv"]= "sans Fe"
            plot_zone_chi2(wavelength, path_to_synth, synth_plot, stardata, spectral_lines,axes=(True,True), size_police=20,size_trace=(1.8, 10),name=name, start=raies.get(str(wavelength))[0],end=raies.get(str(wavelength))[1])
        
    if save:
        with open(save, "w") as fichier:
            json.dump(chi_final, fichier, indent=4, ensure_ascii=False)

# print(raie_Fe.get(str(17012.729)))
# analyse_abu({"22493.67": raie_Fe.get(str(22493.67))}, abu_to_plot=abu_to_plot)



# Afficher une notification avec pync
pync.notify("Votre script a terminé son exécution.", title="Script terminé")