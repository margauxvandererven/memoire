from imports import *
from collections import defaultdict
from scipy import stats

def regression_lineaire(x_vals,y_vals, ax=None):
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

    if ax:
        ax.plot(x_line, y_line, color='firebrick', 
            label=f'Régression: y = ({slope:.3f}±{slope_ci:.3f})x + {intercept:.3f}')

        # Ajout des statistiques sur le graphique
        stats_text = f'R² ajusté = {r_squared_adj:.3f}\n'
        stats_text += f'p-value = {p_value:.3e}\n'
        stats_text += f'Erreur std = {std_residuals:.3f}\n'
        stats_text += f'Normalité p-value = {normality_p_value:.3e}'
        ax.text(0.05, 0.95, stats_text,
                transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        confidence_band = stats.t.ppf(0.975, len(x_vals)-2) * std_residuals
        ax.fill_between(x_line, y_line - confidence_band, y_line + confidence_band, 
            color='gray', alpha=0.1)
        
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

def get_ew_atom(ew_limit, Teff, gamme="IR", particular_element=None):
    if gamme == "IR":
        base_path = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/Linelists/Sophie_IGRINS/"
        files = ["9000-15000_10042024.bsyn",
            "turbo_atoms.20180901_TS2020_transitions_mod_xx_ABO.txt",
            "17000-25000_10042024.bsyn"]
        element_place=11
        wavelen = 14500
    elif gamme == "visible":
        # base_path="/Users/margauxvandererven/MARGAUX/"
        # files = ["nlte_ges_linelist_jmg17feb2022_I_II"]
        base_path="/Users/margauxvandererven/MARGAUX/atoms/"
        files = [
            "3000-4500_Nov11_sansmodif.list",
            "4500-6500_Nov12_GES.list",
            "6500-8900_Nov12_GES.list",
            "8900-10000_Nov12_GES.list"
            ]
        element_place=10
        wavelen = 4000
    data = []
    for file in files:
        with open(base_path + file, "r", encoding="utf8") as f:
            for line in f:
                parts = line.split()
                if len(parts) > 5:
                    element_name = parts[element_place][1].upper() + parts[element_place][2:].lower() + " " + parts[element_place+1]
                    if particular_element and particular_element==element_name:
                        # print(parts)
                        wavelength, excitation_potential, loggf = map(float, parts[:3])
                        ew = 10**(loggf - (5040 / Teff) * excitation_potential)
                        element_name = parts[element_place][1].upper() + parts[element_place][2:].lower() + " " + parts[element_place+1]
                        if ew > ew_limit and wavelength > wavelen:
                            data.append((wavelength, excitation_potential, ew, loggf, element_name))
                    
                    elif particular_element is None:
                        # print(parts)
                        wavelength, excitation_potential, loggf = map(float, parts[:3])
                        ew = 10**(loggf - (5040 / Teff) * excitation_potential)
                        element_name = parts[element_place][1].upper() + parts[element_place][2:].lower() + " " + parts[element_place+1]
                        if ew > ew_limit and wavelength > wavelen:
                            data.append((wavelength, excitation_potential, ew, loggf, element_name))

    return {
        "wavelength": [d[0] for d in data],
        "excitation_potential": [d[1] for d in data],
        "ew": [d[2] for d in data],
        "loggf": [d[3] for d in data],
        "element": [d[4] for d in data],
        "data": data
    }

data = get_ew_atom(ew_limit=1e-9, Teff=4000)["data"]

# lines_BD221742={}
# for d in data:
#      lines_BD221742[d[4]] = []
# for d in data:
#     lines_BD221742.get(d[4]).append(d[0])


# Initialisation du dictionnaire
lines_BD221742 = defaultdict(list)

# Trier les données par d[4] (élément) puis par d[0] (wavelength)
data_sorted = sorted(data, key=lambda x: (x[4], x[0]))

# Dictionnaire temporaire pour stocker la meilleure valeur de d[2] par élément et par plage de d[0]
best_values = {}

for wavelength, excitation_potential, ew, loggf, element in data_sorted:
    # Vérifier s'il existe déjà une valeur proche dans la liste
    if element in best_values:
        close_wavelengths = [w for w in best_values[element] if abs(w - wavelength) < 0.2]

        if close_wavelengths:
            # Vérifier si la nouvelle valeur de d[2] (ew) est plus grande que celle existante
            closest_wavelength = close_wavelengths[0]
            if ew > best_values[element][closest_wavelength]:  
                best_values[element][closest_wavelength] = ew  # Mise à jour avec le plus grand d[2]
        else:
            best_values[element][wavelength] = ew  # Ajouter un nouveau point si aucun proche n'existe
    else:
        best_values[element] = {wavelength: ew}  # Initialiser pour cet élément

# Convertir best_values en dictionnaire {d[4]: [d[0]]}
for element, wavelengths in best_values.items():
    lines_BD221742[element] = list(wavelengths.keys())

# # Affichage du dictionnaire final
# pprint(lines_BD221742)

with open("raies_"+stardata.get("starname")+"_new.txt", "w") as fichier:
    json.dump(lines_BD221742, fichier, indent=4, ensure_ascii=False)


def latex_table(dictionnaire, element):
    sorted_items = sorted(dictionnaire.items())
    latex_table = "\\begin{table}[h!]\n\\centering\n\\small\n\\begin{tabular}{ccc}\n\\hline\n\\hline\n"\
                f"$\\lambda_{{\\mathrm{{{element}}}}}$ & $\\lambda_{{\\mathrm{{min}}}}$ & $\\lambda_{{\\mathrm{{max}}}}$ \\\\\n\\hline\n"

    for key, value in sorted_items:
        latex_table += f"{key:.2f} & {value[0]:.2f}& & {value[1]:.2f}&{value[2]:.2f}& {value[3]:.2f}&{value[4]:.3f}\\\\\n"

    latex_table += "\\end{tabular}\n\\end{table}"

    print(latex_table)

