from astropy.io import fits
from specutils import Spectrum1D
import astropy.units as u
import numpy as np

starname="BD-221742"
spectre="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/H_Sneden/"+starname+"H.cont.fits"

# # Méthode 1 : lecture basique
# hdul = fits.open(spectre)
# data = hdul[0].data  # Accès aux données
# header = hdul[0].header  # Accès aux en-têtes
# hdul.close()

# Méthode 2 : utilisation du context manager (recommandé)
with fits.open(spectre) as hdul:
    hdul.info()
    data = hdul[0].data
    header = hdul[0].header
    for card in header.cards:
        print(card)
    print(data)
    # keywords_of_interest = ['DATE-OBS', 'EXPTIME', 'TELESCOP', 'OBJECT']
    # print("\n=== Mots-clés spécifiques ===")
    # for keyword in keywords_of_interest:
    #     if keyword in hdul[0].header:
    #         print(f"{keyword}: {hdul[0].header[keyword]}")
        
spectrum = Spectrum1D.read(spectre)

flux = spectrum.flux 
wavelength = spectrum.spectral_axis
# print(wavelength)
# print(flux)


with fits.open(spectre) as hdul:
    # Vérifier la structure des données
    print("\n=== Vérification des incertitudes ===")
    
    # 1. Vérifier si un HDU d'erreur existe
    if len(hdul) > 1:
        print("Extensions disponibles:")
        for i, hdu in enumerate(hdul):
            print(f"HDU {i}: {hdu.__class__.__name__}")
            if 'ERROR' in hdu.name.upper() or 'SIGMA' in hdu.name.upper():
                print(f"Trouvé HDU d'erreur: {hdu.name}")
                print("Premiers pixels d'erreur:", hdu.data[:5])
    
    # 2. Vérifier la variance ou erreur dans le header
    header = hdul[0].header
    error_keywords = ['ERROR', 'SIGMA', 'VARIANCE', 'NOISE']
    print("\nMots-clés d'erreur dans le header:")
    for key in header:
        if any(err in key.upper() for err in error_keywords):
            print(f"{key}: {header[key]}")
    
    # 3. Vérifier avec specutils
    spectrum = Spectrum1D.read(spectre)
    if hasattr(spectrum, 'uncertainty'):
        print("\nIncertitudes via Spectrum1D:")
        print("Type:", type(spectrum.uncertainty))
        if spectrum.uncertainty is not None:
            print("Premiers pixels d'incertitude:", spectrum.uncertainty[:5])


def examine_fits(chemin_fichier):
    with fits.open(chemin_fichier) as hdul:
        # 1. Structure générale du fichier
        print("=== STRUCTURE DU FICHIER FITS ===")
        hdul.info()
        
        # 2. Examiner chaque HDU
        for i, hdu in enumerate(hdul):
            print(f"\n=== HDU {i} ===")
            print(f"Type: {hdu.__class__.__name__}")
            
            # Afficher le header
            print("\nHEADER:")
            for card in hdu.header.cards:
                print(f"  {card}")
            
            # Afficher les données
            if hdu.data is not None:
                print("\nDONNÉES:")
                print(f"  Type: {hdu.data.dtype}")
                print(f"  Forme: {hdu.data.shape}")
                print(f"  Premières valeurs: {hdu.data.flatten()[:5]}")
                
                # Si c'est une table
                if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                    print("\nCOLONNES:")
                    for col in hdu.columns:
                        print(f"  - {col.name}: {col.format}")

# Utilisation
spectre = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/spectra/H_Sneden/BD-221742H.cont.fits"
examine_fits(spectre)