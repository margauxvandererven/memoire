from imports import *
from scipy import integrate

def calculate_equivalent_width(wavelength: np.ndarray, normalised_flux: np.ndarray, left_bound: float, right_bound: float) -> float:
    """
    Calculates the equivalent width of a line based on the input parameters
    :param wavelength: Wavelength array
    :param normalised_flux: Normalised flux array
    :param left_bound: Left bound of the line
    :param right_bound: Right bound of the line
    :return: Equivalent width of the line
    """
    # first cut wavelength and flux to the bounds
    try:
        normalised_flux = normalised_flux[np.logical_and(wavelength >= left_bound, wavelength <= right_bound)]
        wavelength = wavelength[np.logical_and(wavelength >= left_bound, wavelength <= right_bound)]
    except TypeError:
        return {"EW": -9999, "error": None}
    
    line_func = interp1d(wavelength, normalised_flux, kind='linear', assume_sorted=True, fill_value=1, bounds_error=False)
    total_area = (right_bound - left_bound) * 1.0   # continuum
    try:
        integration_points = wavelength[np.logical_and.reduce((wavelength > left_bound, wavelength < right_bound))]
        area_under_line, error = integrate.quad(line_func, left_bound, right_bound, points=integration_points, limit=len(integration_points) * 5)
    except ValueError:
        return {"EW": -9999, "error": None}

    return {"EW": total_area - area_under_line, "error": error}

data="../results/Fe_final.txt"
prefix=""

with open(data, "r") as fichier:
    data_lines = json.load(fichier)

ABU=[]
EW=[]
error_EW=[]

for line in data_lines:
    print("Processing line:", line) 
    try : 
        range = str(np.float64(line) - 2) + "-" + str(np.float64(line) + 2)
        abu = data_lines[line][2]
        path_to_synth = "/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/syntspec_TS_local/synth_margaux/"
        synth_name = f"4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_{range}_Feabu_{abu}{prefix}.conv"
        synth_name_blend = f"4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_{range}_Feabu_-20{prefix}.conv"
        synth = syntspec(path_to_synth + synth_name)
        
        synth_blend=syntspec(path_to_synth + synth_name_blend)
        synth_plot={synth_name:str(abu), synth_name_blend:str(-20)}
        
        start=data_lines[line][0]
        end=data_lines[line][1]

        # plot_zone_chi2_simple(np.float64(line), path_to_synth, synth_plot, stardata, lines_BD22, size_police=14,size_trace=(1., 8),name="Fe", start=start,end=end,
        #                         # save="/Users/margauxvandererven/OneDrive  - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/"+str(wavelength)+"/"+str(wavelength)+"_zone", plot=plot
        #                         # save="/Users/margauxvandererven/OneDrive  - Université Libre de Bruxelles/memoire/output/final_Fe/"+str(wavelength)+"/"+str(wavelength)+"_zone_bestfit", plot=True
        #                         )

        ew_line = calculate_equivalent_width(np.array(synth['wavelen']), 
                                   np.array(synth['flux']), 
                                   data_lines[line][0], 
                                   data_lines[line][1])
        ew_blend = calculate_equivalent_width(np.array(synth_blend['wavelen']), 
                                    np.array(synth_blend['flux']), 
                                    data_lines[line][0], 
                                    data_lines[line][1])
        ew = ew_line["EW"] - ew_blend["EW"]        
        
        if ew_line["error"] is not None and ew_blend["error"] is not None:
            error_ew = np.sqrt(ew_line["error"]**2 + ew_blend["error"]**2)
        else:
            error_ew = None        
            
        ew_red = ew / np.float64(line)

        ABU.append(abu)
        EW.append(ew_red)
        error_EW.append(error_ew)
    
    except ValueError as e:
        print(f"Erreur sur la raie {line}: {e}")  # Afficher l'erreur pour debug
        ABU.append(abu)  # On enregistre quand même la raie
        EW.append(None) 
        error_EW.append(None)

for i in EW:
    print(i)


def plot_eq(ABU, ew, chi_final_data, size_police=12, save=None, chi_exc_range=None):
    element_abu="Fe"
    x_vals = []
    y_vals = []
    errors_x=[]

    for i, val in enumerate(ABU):  
        ew_val = np.float64(ew[i]) 
        if ew_val is not None and not np.isnan(ew_val):
            y_vals.append(np.float64(val))  
            x_vals.append(ew_val) 
            errors_x.append(error_EW[i])

    x_vals= np.array(x_vals)
    y_vals=np.array(y_vals)
    errors_x=np.array(errors_x)

    # x_min = 1e-6
    # x_max = 1.75e-5
    # mask = (x_vals > x_min) & (x_vals < x_max)  # Définir x_min et x_max
    # x_vals = x_vals[mask]
    # y_vals = y_vals[mask]
    # errors_x = errors_x[mask]

    errors=np.ones(len(y_vals))*0.1

    print("x_vals", x_vals)
    print("y_vals", y_vals)
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

    ax.errorbar(x_vals, y_vals,yerr=errors,xerr=errors_x, color="gray", label=f"IR : 4000 K", capsize=5,  fmt='o')

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
    
    confidence_band = stats.t.ppf(0.975, len(x_vals)-2) * std_residuals
    ax.fill_between(x_line, y_line - confidence_band, y_line + confidence_band, 
                color='gray', alpha=0.1)
    
    ax.set_xlabel("W$_{\lambda}$/$\lambda$", fontsize=size_police)
    ax.set_ylabel(f"$\\log{{\\epsilon_{{{element_abu}}}}}$", fontsize=size_police)

    ax.xaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    ax.xaxis.set_tick_params(direction='in', length=5, which='minor', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)

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

plot_eq(ABU, EW, data_lines.keys(), size_police=15)