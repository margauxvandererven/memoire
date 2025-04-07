import numpy as np

def syntspec(pathtofile):
    """
    Lecture des fichiers synthétiques
    
    Returns:
        dict: Contient les clés suivantes:
            - 'wavelen': liste des longueurs d'onde
            - 'flux': liste des flux
            - 'header': liste des lignes de commentaires (commençant par ;)
    """
    wavelen = []
    flux = []
    header = []
    with open(pathtofile) as file:
        for line in file:
            line = line.strip()
            if line.startswith(';'):
                header.append(line[1:].strip())  # Enlève le ; et les espaces
            else:
                i = line.split()
                if len(i) >= 2:
                    wavelen.append(float(i[0]))
                    flux.append(float(i[1]))             
    # wavelen = np.array(wavelen)
    # flux = np.array(flux)
    return {
        'wavelen': wavelen, 
        'flux': flux,
        'header': header
    }


def redshift_wavelen(wavelen, v):
    """
    Fonction qui décalent les longueurs d'onde en fonction de redshift
    du à la vitesse radiale de l'étoile
    """
    redshift = []
    c = 299792458 #m/s
    for i in wavelen:
        z = v/c * i
        k = i - z
        redshift.append(k)
    return redshift   


def mediane_5_ancien(liste):
    liste_triee = sorted(liste, reverse=True)
    nb_elements_a_extraire = int(len(liste) * 0.05)
    plus_grands = liste_triee[:nb_elements_a_extraire]
    liste_triee2 = sorted(plus_grands)
    longueur = len(liste_triee2)
    if longueur % 2 != 0:
        mediane_index = longueur // 2
        mediane = liste_triee2[mediane_index]
    else:
        mediane_index_1 = longueur // 2 - 1
        mediane_index_2 = longueur // 2
        mediane = (liste_triee2[mediane_index_1] + liste_triee2[mediane_index_2]) / 2
    
    return mediane

def mediane_5(liste):
    liste_triee = sorted(liste, reverse=True)
    nb_elements_a_extraire = max(1, int(len(liste) * 0.05))

    if nb_elements_a_extraire % 2 == 0:
        nb_elements_a_extraire -= 1

    plus_grands = liste_triee[:nb_elements_a_extraire]
    liste_triee2 = sorted(plus_grands)
    longueur = len(liste_triee2)

    mediane_index = longueur // 2
    mediane = liste_triee2[mediane_index]

    return mediane


def zoom_syntspec(path, filename, wavelength, taille):
    """
    Zoom sur le spectre synthétique pour une longueur d'onde donnée
    """
    wavelenB = syntspec(path+filename)['wavelen']
    fluxB = syntspec(path+filename)['flux']
    h = []
    for j in wavelenB:
        if 0 < j-wavelength:
            h.append(j)
    centre_synt = wavelenB.index(min(h))
    length = int(taille/0.03)
    
    return {"synt_wavelen" : wavelenB[centre_synt-length:centre_synt+length], "synt_flux" : fluxB[centre_synt-length:centre_synt+length]}


def zoom_syntspec2(path, filename, z_wavelen, flux, wavelength):
    """
    Zoom sur le spectre synthétique pour une longueur d'onde donnée, superposé au spectre observé
    """
    wavelenB = syntspec(path+filename)['wavelen']
    fluxB = syntspec(path+filename)['flux']
    k = []
    h = []
    for i in z_wavelen:
        if 0 < i-wavelength:
            k.append(i)
    centre = z_wavelen.index(min(k))

    for j in wavelenB:
        if 0 < j-wavelength:
            h.append(j)
    centre_synt = wavelenB.index(min(h))

    mediane = mediane_5(flux[centre-160:centre+160])

    flux_normalised = []
    for m in flux[centre-160:centre+160] :
        flux_normalised.append(m/mediane)

    return {"z_wavelen" : z_wavelen[centre-160:centre+160], "flux_normalised" : flux_normalised, 
            "synt_wavelen" : wavelenB[centre_synt-2000:centre_synt+2000], "synt_flux" : fluxB[centre_synt-2000:centre_synt+2000]}


def get_nearest(x_query, x_vals):
    return np.abs(x_vals - x_query).argmin()


def normalisation(z_wavelen, flux, wavelength, window_size=20):
    """
    Normalisation du flux observé sur interval de 20 Å. 
    Le flux est divisé par la médiane des 5% plus grandes valeurs sur cet interval.
    """

    window_min = wavelength - window_size/2
    window_max = wavelength + window_size/2

    flux=np.array(flux)
    z_wavelen=np.array(z_wavelen)

    mask = (z_wavelen >= window_min) & (z_wavelen <= window_max)
    window_flux = flux[mask]
    window_wavelen = z_wavelen[mask]

    mediane = mediane_5(window_flux)

    flux_normalised = np.array(window_flux) / mediane
    
    return {"z_wavelen" : window_wavelen, "flux_normalised" : flux_normalised, 
            # "taille" : taille,
            }