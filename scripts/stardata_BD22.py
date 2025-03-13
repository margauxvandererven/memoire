from readmultispec import *
import json
import numpy as np
    
starname = "BD-221742" # nom de l'étoile
starnameb = "BD-221742b" # même étoile, mais normalisation du continu différente 

spectre_visible = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectre_visible/BD221742.coadd"
flux_visible=[]
wavelen_visible=[]

with open(spectre_visible) as f:
    lines=f.readlines()
    for line in lines:
        wavelen_visible.append(np.float64(line.split()[0]))
        flux_visible.append(np.float64(line.split()[1]))

# lecture des spectres fichiers .fits
read_spec_h = readmultispec("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/H_Sneden/"+starname+"H.cont.fits")
read_spec_k = readmultispec("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/K_Sneden/"+starname+"K.cont.fits")
print(read_spec_h['header'])

stardata = {
    "starname"  : starname,
    "wavelen_h" : list(read_spec_h['wavelen']),
    "flux_h" : read_spec_h['flux'], 
    "wavelen_k" : list(read_spec_k['wavelen']), 
    "flux_k" : read_spec_k['flux'], 
    "v_h" : 176.8 * 1e3, # redshift déterminé pour la bande H
    "v_k" : 177.3* 1e3, # redshift déterminé pour la bande K    #
    #  "v_k" : 176.8* 1e3 # redshift déterminé pour la bande K
    "wavelen_visible" : wavelen_visible,
    "flux_visible" : flux_visible
}

read_spec_h_b = readmultispec("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/H_Sneden/"+starnameb+"H.cont.fits")
read_spec_k_b = readmultispec("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/K_Sneden/"+starnameb+"K.cont.fits")

stardatab = {
    "starname"  : starnameb,
    "wavelen_h" : list(read_spec_h_b['wavelen']),
    "flux_h" : read_spec_h_b['flux'], 
    "wavelen_k" : list(read_spec_k_b['wavelen']), 
    "flux_k" : read_spec_k_b['flux'], 
    "v_h" : 176.8 * 1e3, # redshift déterminé pour la bande H
    "v_k" : 177.5* 1e3 # redshift déterminé pour la bande K
}

path = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/syntspec/"+starname+"b/" # chemin vers le dossier des spectres synthétiques

with open("../data_lines/raies_BD-221742_new.txt", "r") as fichier:
    lines_BD22 = json.load(fichier)