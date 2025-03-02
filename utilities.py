from imports import *
from collections import defaultdict


def get_ew_atom(ew_limit, Teff, particular_element=None):
    base_path = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/Linelists/Sophie_IGRINS/"
    files = ["9000-15000_10042024.bsyn",
        "turbo_atoms.20180901_TS2020_transitions_mod_xx_ABO.txt",
        "17000-25000_10042024.bsyn"]
    data = []
    for file in files:
        # print(f"Processing file: {file}")
        with open(base_path + file, "r", encoding="utf8") as f:
            for line in f:
                parts = line.split()
                if len(parts) > 5:
                    element_name = parts[11][1].upper() + parts[11][2:].lower() + " " + parts[12]
                    # Cas 1: particular_element spécifié
                    if particular_element and particular_element==element_name:
                        # print(f"Found line with {atom}: {line.strip()}")
                        wavelength, excitation_potential, loggf = map(float, parts[:3])
                        ew = 10**(loggf - (5040 / Teff) * excitation_potential)
                        element_name = parts[11][1].upper() + parts[11][2:].lower() + " " + parts[12]
                        if ew > ew_limit and wavelength > 14500:
                            data.append((wavelength, excitation_potential, ew, loggf, element_name))
                    
                    # Cas 2: pas de particular_element spécifié
                    elif particular_element is None:
                        wavelength, excitation_potential, loggf = map(float, parts[:3])
                        ew = 10**(loggf - (5040 / Teff) * excitation_potential)
                        element_name = parts[11][1].upper() + parts[11][2:].lower() + " " + parts[12]
                        if ew > ew_limit and wavelength > 14500:
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

