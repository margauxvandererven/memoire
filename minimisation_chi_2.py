import numpy as np
import matplotlib.pyplot as plt
from wavelen_work import *
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import os
import re


def on_click(event):
    if event.xdata is not None and event.ydata is not None:
        coords = f"{event.xdata:.2f}"
        print(coords)
        # Copier les coordonnées dans le presse-papiers
        import pyperclip
        pyperclip.copy(coords)


def get_nearest(x_query, x_vals):
    nearest_index = np.argmin(np.abs(x_vals - x_query))
    return nearest_index


def chi_squared(synth, observed):
    return np.sum((synth - observed) ** 2 / observed)

def zone_chi2(path, synthetics, stardata, k, spectral_lines,chi_final, name=None, start=None, end=None, plot=None, save=None):
    if k < 18500:
        j="h"
    else :
        j="k"
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + j), stardata.get("v_" + j)), stardata.get("flux_" + j), k)
    wavelength_index = get_nearest(k, np.array(normal["z_wavelen"]))
    tolerance = 0.05
    index_closest_with_tolerance = next((i for i, x in enumerate(normal['flux_normalised'][wavelength_index:]) if abs(x - 1) <= tolerance), None)
    # Longueurs d'onde observées
    observed_wavelengths = np.array(normal['z_wavelen'][wavelength_index - index_closest_with_tolerance +1 : wavelength_index + index_closest_with_tolerance ])
    # print(observed_wavelengths[0], observed_wavelengths[-1])
    # Flux observé : portion autour de l'indice central
    observed = np.array(normal['flux_normalised'][wavelength_index - index_closest_with_tolerance+1:wavelength_index + index_closest_with_tolerance ])
    
    if start is not None and end is not None:
        observed_wavelengths = np.array(normal['z_wavelen'][get_nearest(start,np.array(normal['z_wavelen'])):get_nearest(end,np.array(normal['z_wavelen']))])
        observed = np.array(normal['flux_normalised'][get_nearest(start,np.array(normal['z_wavelen'])):get_nearest(end,np.array(normal['z_wavelen']))])

    chi_final[k]=[observed_wavelengths[0], observed_wavelengths[-1]]
    
    taille = (observed_wavelengths[-1]-observed_wavelengths[0])/2
    # Création de la liste pour stocker les spectres 
    synthetic_spectra = []
    observed_spectra =[]
    observed_spectra_w =[]
    synthetic_spectra_w = []

#tracé du graphe
    if plot is True:
        f = plt.figure(figsize=(8,4))
        gs = f.add_gridspec(1, hspace=0.2)
        ax = gs.subplots(sharex=False, sharey=True)

        ax.set_xlim(k-3,k+3)
        ax.set_ylim(0., 1.4)
        ax.set_xlabel("Longueur d'onde (Å)", fontsize  = 10)

        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))

        # Récupération de toutes les longueurs d'onde dans une liste unique avec les éléments associés
        all_wavelengths = []
        for element, wavelengths in spectral_lines.items():
            for wavelength in wavelengths:
                all_wavelengths.append((wavelength, element))

        # Tri de la liste des longueurs d'onde
        all_wavelengths.sort()  # Trie par longueur d'onde en ordre croissant

        # Boucle pour tracer les lignes et afficher les éléments
        text_height_base = 1.2
        text_height = text_height_base  # Initialisation de la hauteur actuelle

        for i, (z, element) in enumerate(all_wavelengths):
            # Détection de la proximité avec la longueur d'onde précédente
            if i > 0 and abs(z - all_wavelengths[i - 1][0]) < 0.7:
                text_height += 0.07  # Si trop proche du précédent, augmenter la hauteur de 0.05
            else:
                text_height = text_height_base  # Réinitialisation de la hauteur si suffisamment éloigné

            # Vérification des conditions pour tracer les lignes
            if k - 11 <= z <= k + 11:
                # Tracé de la ligne verticale
                ax.axvline(x=z, ymin=0.7, ymax=0.75, color='black', linewidth=0.5)
                # Ajout du texte de l'élément avec hauteur adaptée
                ax.text(z, text_height, s=element, color='black', fontsize=10, ha='center')
            
            if k - 3 <= z <= k + 3:
                # Tracé de la ligne verticale
                ax.axvline(x=z,   ymin=0.7, ymax=0.75, color='black', linewidth=0.5)
                # Ajout du texte de l'élément avec hauteur adaptée
                ax.text(z, text_height, s=element, color='black', fontsize=10, ha='center')
        
        ax.scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=5,facecolors='none', color='black', label="Spectre observé : " + stardata.get("starname"))
        
        if name is not None:
            ax.axvline(x=k,   ymin=0.7, ymax=0.75, color='black', linewidth=0.5)
            ax.text(k, 1.2, s=name, color='black', fontsize=10, ha='center')
        
        for synth in synthetics : 
            AX2 = syntspec(path+synth)
            ax.plot(AX2['wavelen'], AX2['flux'], linewidth = 1, label="Synth : "+synthetics.get(synth))

        # ax.plot(observed_spectra_w[0],observed_spectra[0])

        ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
        ax.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
        ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
        
        ax.tick_params(axis = 'both', labelsize = 10)
        ax.legend(loc = "lower left", fontsize = 10)
        ax.axvspan(observed_wavelengths[0], observed_wavelengths[-1], color='bisque', alpha=0.3)
        ax.set_ylabel("Flux normalisé", fontsize  = 10)
        plt.show()

        if save is not None:
            plt.savefig(save, dpi=400)


def chi_2(path, synthetics, stardata, k, spectral_lines,chi_final,start, end, name=None, plot=None, save=None, size_police=None, axes=None):
    if k < 18500:
        j="h"
    else :
        j="k"
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + j), stardata.get("v_" + j)), stardata.get("flux_" + j), k)
    wavelength_observed = np.array(normal['z_wavelen'])
    flux_observed = np.array(normal['flux_normalised'])

    mask_observed = (wavelength_observed >= start) & (wavelength_observed <= end)
    wavelength_observed_filtered = wavelength_observed[mask_observed]
    flux_observed_filtered = flux_observed[mask_observed]

    chi_final[k]=[start, end]
    
    synthetic_spectra = []
    observed_spectra_interpolated =[]
    observed_spectra_w_interpolated =[]
    synthetic_spectra_w = []
    chi_squared_values = []

    for syntha in synthetics:
            synt_flux = np.array(zoom_syntspec(path, syntha, k, 10)["synt_flux"])
            synt_w = np.array(zoom_syntspec(path, syntha, k, 10)["synt_wavelen"])

            mask_synthetic = (synt_w >= start) & (synt_w <= end)

            wavelength_synthetic_filtered = synt_w[mask_synthetic]
            flux_synthetic_filtered = synt_flux[mask_synthetic]

            synthetic_spectra.append(flux_synthetic_filtered)
            synthetic_spectra_w.append(wavelength_synthetic_filtered)

            interpolator = interp1d(wavelength_observed_filtered, flux_observed_filtered, kind='linear', fill_value="extrapolate")
            flux_observed_interpolated = interpolator(wavelength_synthetic_filtered)

            chi2 = chi_squared(flux_synthetic_filtered, flux_observed_interpolated)
            chi_squared_values.append(chi2)

            observed_spectra_interpolated.append(flux_observed_interpolated)
            observed_spectra_w_interpolated.append(wavelength_synthetic_filtered)


    items = list(synthetics.items())

    if plot is True:
        f = plt.figure(figsize=(12,8))
        gs = f.add_gridspec(2, hspace=0.2)
        ax = gs.subplots(sharex=False, sharey=True)
        cid = f.canvas.mpl_connect('button_press_event', on_click)

        ax[0].set_xlim(k-10,k+10)
        ax[1].set_xlim(k-2,k+2)
        ax[0].set_ylim(0.0, 1.4)
        ax[1].set_ylim(0.0, 1.4)
        if size_police is None:
            size_police = 10

        if axes is True:
            ax[0].set_ylabel("Flux normalisé", fontsize  = size_police)
            ax[1].set_xlabel("Longueur d'onde (Å)", fontsize  = size_police)

        ax[0].xaxis.set_major_locator(MultipleLocator(5))
        ax[0].xaxis.set_minor_locator(MultipleLocator(1))
        ax[1].xaxis.set_major_locator(MultipleLocator(1))
        ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))

        all_wavelengths = []
        for element, wavelengths in spectral_lines.items():
            for wavelength in wavelengths:
                all_wavelengths.append((wavelength, element))

        all_wavelengths.sort()  

        text_height_base = 1.2
        text_height = text_height_base  

        for i, (z, element) in enumerate(all_wavelengths):
            if i > 0 and abs(z - all_wavelengths[i - 1][0]) < 0.7:
                text_height += 0.07  
            else:
                text_height = text_height_base  

            if k - 11 <= z <= k + 11:
                ax[0].axvline(x=z, ymin=0.7, ymax=0.75, color='black', linewidth=1)
                ax[0].text(z, text_height, s=element, color='black', fontsize=size_police, ha='center')
            
            if k - 3 <= z <= k + 3:
                ax[1].axvline(x=z,   ymin=0.7, ymax=0.75, color='black', linewidth=1)
                ax[1].text(z, text_height, s=element, color='black', fontsize=size_police, ha='center')
        
        for ax_ in ax :
            ax_.scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=8,facecolors='none', color='black', 
                    #    label="Spectre observé : " + stardata.get("starname")
                       )
            if name is not None:
                ax_.axvline(x=k,   ymin=0.7, ymax=0.75, color='black', linewidth=0.5)
                ax_.text(k, 1.2, s=name, color='black', fontsize=size_police, ha='center')
            
            for synth in synthetics : 
                AX2 = syntspec(path+synth)
                ax_.plot(AX2['wavelen'], AX2['flux'], linewidth = 1.5, label=synthetics.get(synth))
            ax_.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
            ax_.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
            ax_.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
            ax_.xaxis.get_major_formatter().set_scientific(False)
            ax_.xaxis.get_major_formatter().set_useOffset(False)
            ax_.tick_params(axis = 'both', labelsize = size_police)
            ax_.axvspan(start, end, color='bisque', alpha=0.2)
        ax[0].legend(fontsize = size_police, framealpha=0.8, facecolor='whitesmoke', markerscale=0.2,edgecolor='whitesmoke',numpoints=5)
        plt.tight_layout()
        plt.show()
        if save is not None:
            plt.savefig(save + '.png', dpi=400, transparent=True)
            # plt.savefig(save + '.pgf', backend='pgf')

    return {"chi_squared_values":chi_squared_values}


def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

def chi_plot_ABU(abu, chi_squared,element, k, raie, raie_OH):
    # plt.switch_backend('pgf')
    abu = np.array(abu)
    # Perform quadratic fitting
    params, _ = curve_fit(quadratic, abu, chi_squared)
    a, b, c = params  # coefficients of the quadratic fit
    fit_x = np.linspace(abu.min(), abu.max(), 100)
    fit_y = quadratic(fit_x, a, b, c)

    # Find minimum of the quadratic fit
    min_log_e = -b / (2 * a)  # Vertex of the parabola
    min_chi_squared = quadratic(min_log_e, a, b, c)

    # Plot
    f = plt.figure(figsize=(8, 6))
    gs = f.add_gridspec(1, hspace=0.2)
    ax = gs.subplots(sharex=False, sharey=True)
    ax.scatter(abu, chi_squared, color="darkblue", marker="x", label=f"Raie de {raie} en {k} Å")  # Data points
    ax.plot(fit_x, fit_y, color="lightgray")  # Quadratic fit line
    ax.scatter(min_log_e, min_chi_squared, color="red", marker="x", label=f"Minimum en {min_log_e:.2f}")  # Minimum point
    ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
    ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
    ax.tick_params(axis = 'both', labelsize = 16)
    # Labeling
    ax.set_xlabel(f"$\\log \\epsilon_{{\\mathrm{{{element}}}}}$", fontsize=16)
    ax.set_ylabel(r"$\chi^2$", fontsize=16)
    plt.legend(loc="upper left", fontsize=16)
    plt.show()
    # plt.savefig('graphique_test.pgf')
    raie_OH[k] = min_log_e

def interpolated_chi_squared(x, abundances, macroturbulences, chi_squared_values):
    ab, mt = x
    return griddata((abundances, macroturbulences), chi_squared_values, (ab, mt), method='cubic')

# Fonction de minimisation et visualisation
def double_chi(abundances, macroturbulences, chi_squared_values):
    ab_min, ab_max = min(abundances), max(abundances)
    mt_min, mt_max = min(macroturbulences), max(macroturbulences)

    # Générer une grille fine pour la visualisation
    grid_ab, grid_mt = np.mgrid[ab_min:ab_max:100j, mt_min:mt_max:100j]
    grid_chi_squared = griddata((abundances, macroturbulences), chi_squared_values, (grid_ab, grid_mt), method='cubic')

    # Estimation initiale pour la minimisation
    initial_guess = [np.mean(abundances), np.mean(macroturbulences)]
    
    # Minimisation avec interpolation du chi carré
    result = minimize(
        interpolated_chi_squared, initial_guess, 
        args=(abundances, macroturbulences, chi_squared_values),  # Arguments supplémentaires
        bounds=[(ab_min, ab_max), (mt_min, mt_max)]
    )

    # Récupérer les valeurs minimales estimées
    best_abundance = result.x[0]
    best_macroturbulence = result.x[1]
    min_chi_squared_value = result.fun

    print("Meilleure abondance estimée :", best_abundance)
    print("Meilleure macroturbulence estimée :", best_macroturbulence)
    print("Valeur minimale de χ² estimée :", min_chi_squared_value)

    # Visualisation
    plt.figure(figsize=(10, 8))
    plt.contourf(grid_ab, grid_mt, grid_chi_squared, levels=20, cmap='magma')
    plt.colorbar(label='$\chi^2$')
    plt.scatter(abundances, macroturbulences, c=chi_squared_values, cmap='magma', edgecolor='orange', label='Données')
    plt.plot(best_abundance, best_macroturbulence, 'rx', markersize=10, label='Min $\chi^2$ interpolé')
    plt.xlabel("Abondance")
    plt.ylabel("Macroturbulence")
    plt.legend()
    plt.show()


def plot_chi_squared_MAC(path, stardata,lines_BD22, raie_propre, dossier_element, chi_final):
    """
    path : chemin du dossier vers spectres synthétiques
    stardata : données de l'étoile
    lines_BD22 : raies atomiques de tout le spectre
    raie_propre : raie étudiée
    dossier_element : dossier vers spectres synthétiques de la raie étudiée
    """
    abondance = [] # on crée liste pour stocker les abondances suite à la lecture des fichiers synthétiques
    macroturbulence = []  # on crée liste pour stocker les macroturb suite à la lecture des fichiers synthétiques
    synth = {}
    target_path = path+dossier_element+"/"
    for filename in os.listdir(target_path):
        # Vérifier si c'est un fichier (et non un dossier) et ne pas inclure .DS_Store
        if os.path.isfile(os.path.join(target_path, filename)) and filename != ".DS_Store":
            match = re.search(r'MAC_-(.*?)_', filename) 
            match2 = re.search(r'abu_(.*?).conv', filename)          
            if match and match2:
                chars_after_MAC = match.group(1)
                chars_after_ABU = match2.group(1)

            synth[filename] = f"log $\\epsilon_{{\\mathrm{{{dossier_element}\\; I}}}}$ = {chars_after_ABU} & $v_{{\\mathrm{{macro}}}}$ = {chars_after_MAC} $km s^{{-1}}$"
            abondance.append(float(chars_after_ABU))
            macroturbulence.append(float(chars_after_MAC))

    # Calculer les valeurs de chi_squared pour chaque raie
    chi_squared_values = None  # Initialiser au cas où aucune valeur ne soit calculée
    for raie in raie_propre:
        for k in raie_propre[raie]:
            result = chi_2(target_path, synth, stardata, k, lines_BD22,chi_final, name=raie)
            chi_squared_values = result.get("chi_squared_values", [])  # Assurez-vous de récupérer chi_squared_values

    # Appeler double_chi si chi_squared_values a bien été défini
    if chi_squared_values is not None:
        double_chi(abondance, macroturbulence, chi_squared_values)
    else:
        print("Erreur: chi_squared_values n'a pas été calculé.")



def plot_chi2_simple_MAC(path, dossier_element, stardata, raie_propre, lines_BD22,chi_final, target_pairs=None):
    """
    Affiche les graphes pour les spectres synthétiques qui correspondent aux paires d'abondance et de macroturbulence cibles.

    target_pairs : Liste de tuples avec des paires d'(abondance, macroturbulence). Exemple : [(6.5, 10.5), (7.0, 9.0)]
    """
    synth = {}
    target_path = path+dossier_element+"/"
    
    for filename in os.listdir(target_path):
        # Vérifier si c'est un fichier (et non un dossier) et ne pas inclure .DS_Store
        if os.path.isfile(os.path.join(target_path, filename)) and filename != ".DS_Store":
            # Rechercher les valeurs de MAC et d'abondance dans le nom du fichier
            match = re.search(r'MAC_-(.*?)_', filename)
            
            match2 = re.search(r'abu_(.*?).conv', filename)
            
            if match and match2:
                chars_after_MAC = match.group(1)
                chars_after_ABU = match2.group(1)

                # Vérifier si la paire (abondance, macroturbulence) est dans target_pairs
                if target_pairs is None or (float(chars_after_ABU), float(chars_after_MAC)) in target_pairs:
                    synth[filename] = f"log $\\epsilon_{{\\mathrm{{{dossier_element}\\; I}}}}$ = {chars_after_ABU} & $v_{{\\mathrm{{macro}}}}$ = {chars_after_MAC} $km s^{{-1}}$"
                    
    # Parcourir les raies et appeler la fonction chi_2 pour chaque raie spécifique
    for raie in raie_propre:
        for k in raie_propre[raie]:
            chi_2(target_path, synth, stardata, k, lines_BD22,chi_final, plot=True)


def plot_chi2_simple_ABU(path, dossier_element, stardata, raie_propre, lines_BD22,start=None,end=None, target_pairs=None):

    """

    """
    synth = {}
    target_path = path+dossier_element+"/"
    
    for filename in os.listdir(target_path):
        # Vérifier si c'est un fichier (et non un dossier) et ne pas inclure .DS_Store
        if os.path.isfile(os.path.join(target_path, filename)) and filename != ".DS_Store":
            # Rechercher les valeurs de MAC et d'abondance dans le nom du fichier
            
            match = re.search(r'abu_(.*?).conv', filename)
            
            if match:
                chars_after_ABU = match.group(1)

                # Vérifier si la paire (abondance, macroturbulence) est dans target_pairs
                if target_pairs is None or (float(chars_after_ABU)) in target_pairs:
                    synth[filename] = f"log $\\epsilon_{{\\mathrm{{{dossier_element}\\; I}}}}$ = {chars_after_ABU}"
                    
    # Parcourir les raies et appeler la fonction chi_2 pour chaque raie spécifique
    for raie in raie_propre:
        for k in raie_propre[raie]:
            chi_2(target_path, synth, stardata, k, lines_BD22,start,end, plot=True)

def abu_plot(raies_element, element_abu, save=None):
    x_vals = list(raies_element.keys())
    y_vals = list(raies_element.values())

    f = plt.figure(figsize=(10, 5))
    gs = f.add_gridspec(1)
    ax = gs.subplots(sharex=False, sharey=True)
    
    # Scatter plot for raies
    ax.scatter(x_vals, y_vals, color="darkblue", label=f"$\\log_{{\\epsilon_{{{element_abu}}}}}$", marker="x")
    ax.set_xlabel("$\\lambda$ (Å)", fontsize=12)
    ax.set_ylabel(f"$\\log_{{\\epsilon_{{{element_abu}}}}}$", fontsize=12)
    
    # Set tick parameters
    ax.xaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    ax.xaxis.set_tick_params(direction='in', length=5, which='minor', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    
    # Set limits
    # ax.set_xlim(14700, 18000)
    # ax.set_ylim(8.5, 8.7)

    # Calculate mean and standard deviation, rounded to two decimal places
    y_mean = round(np.mean(y_vals), 2)
    y_std = round(np.std(y_vals), 2)

    # Filter points within the band
    y_min_band = y_mean - y_std
    y_max_band = y_mean + y_std
    filtered_y_vals = [y for y in y_vals if y_min_band <= y <= y_max_band]
    filtered_x_vals = [x_vals[i] for i, y in enumerate(y_vals) if y_min_band <= y <= y_max_band]

    # Recalculate the mean and standard deviation for filtered points
    filtered_mean = round(np.mean(filtered_y_vals), 2)
    filtered_std = round(np.std(filtered_y_vals), 2)

    # Plot horizontal line for the filtered mean
    ax.hlines(y=y_mean, xmin=14500, xmax=23500, color="red", linestyle="--", linewidth=1, label=f"Moyenne = {y_mean:.2f}")

    # Plot original shaded error band
    ax.fill_between(
        [14500, 23500], 
        y_mean - y_std, 
        y_mean + y_std, 
        color="gray", 
        alpha=0.1, 
        label=f"$\sigma$ = ±{y_std:.2f}"
    )

    # Display legend
    plt.legend()
    plt.show()
    if save is not None:
        plt.savefig(save, dpi=400)

    # Return filtered mean and standard deviation
    return filtered_mean, filtered_std