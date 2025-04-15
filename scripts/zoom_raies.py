import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from wavelen_work import *
import os 
from scipy.signal import find_peaks
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# lines = {
#     "S": [15478.48, 22507.60],
#     "Na": [16373.87 , 16388.85 , 22056.43 , 22083.66  , 23379.14],
#     "Al": [16718.97 , 16750.60 , 16763.37 , 17699.05 , 21093.08 , 21163.80 , 21208.18], 
#     "Zn": [16505.18],
#     "Si": [16434.93, 20804.20 , 20890.37  , 20926.14],
#     "K": [15163.09 , 15168.40],
#     "Ca": [16150.76 , 16155.24 , 16157.36 ,  20962.57 , 20972.53 , 21113.90 , 22607.94 , 22624.96 , 22626.72 , 22651.18 , 22653.58],
#     "Sc": [21730.36 , 21812.24 , 21842.70],
#     "Ti": [16330.54 , 21149.62],
#     "V": [15924.81],
#     "Cr": [15680.06 , 15860.21 , 17708.73 ], 
#     "Mn": [15217.74 , 15262.49],
#     "Co": [16757.64],
#     "Ni": [16310.50 , 16363.09 , 16815.46 , 16818.74 , 16867.28 , 17306.52 , 20957.14 , 21167.93 , 21570.06 , 21945.50 , 22596.93 ],
#     "Cu": [16005.64, 16638.98], 
#     "Zn": [16505.18],
#     "Y": [21260.45 , 22543.84],
#     "Ce": [15277.65 , 15829.83 , 15977.12 , 16595.18],
#     "Nd": [15368.14 , 16053.63 , 16262.04],
#     "Yb": [16498.40],
#     "F": [22699.49 , 22714.59 ,  22778.25 , 22826.86 , 22886.73 , 22957.94],
#     "Mg" : [21059.76, 21060.89, 21458.87]
# }

def find_peaks_element(path, fichier1, fichier2, plot=None):
    fluxA = np.array(syntspec(path+fichier1)['flux'])
    fluxB = np.array(syntspec(path+fichier2)['flux'])
    wavelength = syntspec(path+fichier1)['wavelen']

    residu = fluxA - fluxB
    # print(residu)

    ta = []
    for i in range(len(residu)):
        if np.abs(residu[i]) > 0.01:
            ta.append(wavelength[i])

    # print(len(ta))

    longueur_onde = np.array(syntspec(path+fichier1)['wavelen'])
    flux = 1-np.abs(residu)

    minima, _ = find_peaks(-flux, prominence=0.05) #seul les pics qui dépassent leur voisins de plus de 5% sont pris en compte
    if plot:
        plt.figure(figsize=(8, 5))
        plt.plot(longueur_onde, flux, label="Spectre")
        plt.plot(longueur_onde[minima], flux[minima], "ro", label="Pics d'absorption")
        plt.xlabel("Longueur d'onde (nm)")
        plt.ylabel("Flux")
        plt.title("Détection des pics d'absorption")
        plt.legend()
        plt.show()

    return {"minima" : longueur_onde[minima]}



def validate_wavelengths(wavelengths, path, synthetics, stardata,taille_zoom, spectral_lines):
    """
    Permet de valider interactivement une liste de longueurs d'onde.
    
    Returns:
        list: Liste des longueurs d'onde validées
    """
    validated_wavelengths = []
    
    for k in wavelengths:
        is_valid = zoom_lines({"":[k]},path, synthetics, stardata, taille_zoom, spectral_lines, interactive=True
        )
        
        if is_valid:
            validated_wavelengths.append(k)
            print(f"Raie {k} Å ajoutée à la liste")
        else:
            print(f"Raie {k} Å ignorée")
            
    return validated_wavelengths


def on_click(event):
    if event.xdata is not None and event.ydata is not None:
        coords = f"{event.xdata:.2f}"
        print(coords)
        # if wavelength_range:
        #     wavelength_range.append(coords)
        # Copier les coordonnées dans le presse-papiers
        import pyperclip
        pyperclip.copy(coords)




def plot_lines_interactive(path, synthetics, stardata, k, taille, spectral_lines, gamme="IR"):
    """
    Version interactive de plot_lines utilisant plotly
    """
    # Préparation des données
    if gamme == "IR":
        normal = normalisation(redshift_wavelen(stardata.get("wavelen_h"), 
                             stardata.get("v_h")), stardata.get("flux_h"), k)
    elif gamme == "visible":
        normal = normalisation(stardata.get("wavelen_visible"), 
                             stardata.get("flux_visible"), k, window_size=100)
    
    # Création de la figure interactive
    fig = go.Figure()

    # Ajout du spectre observé
    fig.add_trace(go.Scatter(
        x=normal['z_wavelen'],
        y=normal['flux_normalised'],
        mode='lines',
        name=f"Spectre observé : {stardata.get('starname')}",
        line=dict(color='black', dash='dash')
    ))

    # Ajout des spectres synthétiques
    for mol in synthetics:
        AXb = syntspec(path + mol)
        fig.add_trace(go.Scatter(
            x=AXb['wavelen'],
            y=AXb['flux'],
            mode='lines',
            name=f"Synth: {synthetics.get(mol)}"
        ))

    # Configuration de la mise en page
    fig.update_layout(
        title=f"Spectre centré sur {k} Å",
        xaxis_title="Longueur d'onde (Å)",
        yaxis_title="Flux normalisé",
        xaxis=dict(
            range=[k - taille, k + taille],
            tickmode='linear',
            dtick=1
        ),
        yaxis=dict(
            range=[-0.2 if gamme == "visible" else 0, 1.3],
        ),
        hovermode='x',
        showlegend=True
    )

    # Ajout des lignes verticales pour les éléments
    for element, wavelengths in spectral_lines.items():
        for z in wavelengths:
            if k - taille <= z <= k + taille:
                fig.add_vline(
                    x=z,
                    line_width=1,
                    line_dash="dash",
                    annotation_text=element,
                    annotation_position="top"
                )

    return fig

# Utilisation:
# fig = plot_lines_interactive(path, synthetics, stardata, k, taille_zoom, spectral_lines)
# fig.show()


def plot_lines(ax,f, l, path, synthetics, stardata, n, k, i, m, j, taille, spectral_lines, save,gamme, range=None, plot=True, mol_band=None):
    """
    Ouvre et trace le spectre observé normalisé. Trace les spectres synthétiques "synthetics". Le tout centré sur les raies séléctionnées
    avec une largeur demandée (taille)
    ax : le plot
    l : ?
    path : chemin des spectres synthétiques
    synthetics : spectre synthétiques
    stardata : données de l'étoile
    n : nombre d'éléments dans la liste de raies
    k : raie de centrage
    i : nom de la raie
    m : bande h ou k
    j : "", "2" ou "3"
    taille : taille du zoom
    """
    colors = ['#3e92f7', 
            #   '#faaa2e',
               '#db2f2f','#a39cf3', '#92c124']
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)

    fontsize=14
    linesize=1.0

    if gamme == "IR":
        normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + m), stardata.get("v_"+ m)), stardata.get("flux_" + m), k)
    elif gamme == "visible":
        normal = normalisation(stardata.get("wavelen_visible"), stardata.get("flux_visible"), k, window_size=100)
       
    if n == 1:
        if gamme == "IR":
            ax.scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=5,facecolors='none', color='black', label="Spectre observé : " + stardata.get("starname"))
        if gamme == "visible":
            ax.plot(normal['z_wavelen'], normal['flux_normalised'], color='black', label="Spectre observé : " + stardata.get("starname"), linewidth=0.6, linestyle='--', alpha=0.8, dashes=(10, 4)) 
        for mol in synthetics:
            AXb = syntspec(path + mol + j)
            ax.plot(AXb['wavelen'], AXb['flux'], label= "Synth: " + str(synthetics.get(mol)), linewidth=1.)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_tick_params(direction='in', length=10, which='major')
        ax.xaxis.set_tick_params(direction='in', length=6, which='minor')
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.get_major_formatter().set_scientific(False)
        ax.xaxis.get_major_formatter().set_useOffset(False)

        all_wavelengths = []
        for element, wavelengths in spectral_lines.items():
            for wavelength in wavelengths:
                all_wavelengths.append((wavelength, element))
        # Tri de la liste des longueurs d'onde
        all_wavelengths.sort()
        text_height_base = 1.2
        text_height = text_height_base  # Initialisation de la hauteur actuelle

        for i, (z, element) in enumerate(all_wavelengths):
            # Détection de la proximité avec la longueur d'onde précédente
            if i > 0 and abs(z - all_wavelengths[i - 1][0]) < 0.7:
                text_height += 0.07  # Si trop proche du précédent, augmenter la hauteur de 0.05
            else:
                text_height = text_height_base  # Réinitialisation de la hauteur si suffisamment éloigné

            # Vérification des conditions pour tracer les lignes
            if k - taille <= z <= k + taille:
                # Tracé de la ligne verticale
                ax.axvline(x=z, ymin=0.8, ymax=0.85, color='black', linewidth=0.5)
                # Ajout du texte de l'élément avec hauteur adaptée
                ax.text(z, text_height, s=element, color='black', fontsize=10, ha='center')

        # ax.text(k, 1.2, s=i, color='black', fontsize=12, ha='center')
        # ax.axvline(x=k, ymin=0.8, ymax=0.9, color='black', linewidth=0.5)
        ax.set_xlim(k - taille, k + taille)
        ax.set_xlim(k - taille, k + taille)
        if gamme == "IR":
            ax.set_ylim(0., 1.3)
        elif gamme == "visible":
            ax.set_ylim(-0.4, 1.1)
        # ax.axhline(y=0.3,xmin= k - 10,xmax= k+10, color='red', linewidth=0.5)
        ax.legend(loc='lower right')
        if range is not None:
            def on_click_range(event):
                if event.xdata is not None:
                    wavelength = round(event.xdata, 3)
                    if k not in range:
                        range[k] = [wavelength]
                    elif len(range[k]) < 2:
                        range[k].append(wavelength)
                        range[k].sort()
                        print(f"Plage pour {k} Å : {range[k]}")

            print(f"\nSélection de la plage pour la raie {k} Å")
            print("Cliquez deux fois pour définir le début et la fin de la plage")
                        
            f.canvas.mpl_connect('button_press_event', on_click_range)
        else:
            cid = f.canvas.mpl_connect('button_press_event', on_click)
        if save is not None:
            # Extraire le chemin du dossier depuis le chemin complet
            save_dir = os.path.dirname(save)
            
            # Créer le dossier s'il n'existe pas
            if save_dir and not os.path.exists(save_dir):
                os.makedirs(save_dir)
                
            plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)
        if plot:
            plt.show()
        else:
            plt.close()
    else:
        print("bonjour")
        
        if gamme == "IR":
            ax[l].scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=5,facecolors='none', color='black', label="Spectre observé : " + stardata.get("starname"))
        if gamme == "visible":
            ax[l].plot(normal['z_wavelen'], normal['flux_normalised'], color='black', label="Spectre observé : " + stardata.get("starname") , linewidth=linesize)
        
        for idx, mol in enumerate(synthetics):
            print("rebonjour")
            AXb = syntspec(path + mol + j)
            ax[l].plot(AXb['wavelen'], AXb['flux'], label="Synth: " + str(synthetics.get(mol)), linewidth=linesize,color=colors[idx % len(colors)] )
        
        ax[l].hlines(y=1,xmin= k - taille,xmax= k+taille, color='gray', linewidth=0.5, alpha=0.8)
        ax[l].xaxis.set_major_locator(MultipleLocator(5))
        ax[l].xaxis.set_minor_locator(MultipleLocator(1))
        ax[l].xaxis.set_tick_params(direction='in', length=10, which='major')
        ax[l].xaxis.set_tick_params(direction='in', length=6, which='minor')
        ax[l].yaxis.set_tick_params(direction='in', length=10, which='major')
        ax[l].yaxis.set_tick_params(direction='in', length=6, which='minor')
        ax[l].xaxis.set_major_formatter(ScalarFormatter())
        ax[l].xaxis.get_major_formatter().set_scientific(False)
        ax[l].xaxis.get_major_formatter().set_useOffset(False)
        ax[l].tick_params(axis = 'both', labelsize = fontsize)
        for element, wavelengths in spectral_lines.items():
            for z in wavelengths:
                if k - taille <= z <=  k + taille:
                    # Tracé de la ligne verticale
                    ax[l].axvline(x=z, ymin=0.8, ymax=0.9, color='darkgray', linewidth=1.)
                    # Ajout du texte de l'élément au-dessus de la ligne
                    ax[l].text(z, 1.2, s=element, color='black', fontsize=fontsize, ha='center')
        ax[l].axvline(x=k, ymin=0.8, ymax=0.9, color='black', linewidth=0.5)
                    # Ajout du texte de l'élément au-dessus de la ligne
        ax[l].text(k, 1.1, s=i, color='black', fontsize=fontsize, ha='center')
        ax[l].set_xlim(k - taille, k + taille)

        if mol_band is not None:
            for band in mol_band:
                if band[0] < k + taille and band[1] > k - taille:
                    ax[l].axvspan(band[0], band[1], color='lightgray', alpha=0.7)
                    ax[l].text((band[0] + band[1]) / 2, 1.05, "CN", color='black', fontsize=fontsize, ha='center')

        if gamme == "IR":
            ax[l].set_ylim(0., 1.3)
        elif gamme == "visible":
            ax[l].set_ylim(-0.2, 1.2)

        if l == n-1:
            ax[l].legend(loc='lower right', fontsize=fontsize)
            ax[l].set_xlabel("Longueur d'onde [Å]", fontsize=fontsize)
        if l == n//2: 
            ax[l].set_ylabel("Flux normalisé", fontsize=fontsize)
        
        


        cid = f.canvas.mpl_connect('button_press_event', on_click)
        
        # if save is not None:
        #     save_dir = os.path.dirname(save)
        #     if save_dir and not os.path.exists(save_dir):
        #         os.makedirs(save_dir)
        #     plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)

        # if plot:
        #     plt.show()
        # else:
        #     plt.close()
         

def zoom_lines(lines, path, synthetics, stardata, taille_zoom, spectral_lines, gamme="IR", save=None, interactive=False, raies_validees=None, range=None, plot=True, mol_band=None):
    """
    lines : raies que l'on souhaite observer
    path : chemin vers spectres synthétiques
    synthetics : spectres synthétiqyes
    stardata : données de l'étoile
    taille_zoom : taille du zoom sur chaque raie
    """
    # raies_validees = []
    for i in list(lines.keys()):
        n = len(lines.get(i))
        y_size = n*3
        f = plt.figure(figsize=(20,y_size))
        gs = f.add_gridspec(n, hspace=0.2)
        if n == 1:
            ax = gs.subplots()  # Sans sharex et sharey pour un seul subplot
            ax = [ax]  # Convertir en liste pour uniformiser le traitement
        else:
            ax = gs.subplots(sharex=False, sharey=True)
        for l, k in enumerate(lines.get(i)):
            print(l)
            # if gamme == "IR":
            if k < 18500:
                plot_lines(ax, f, l, path, synthetics, stardata, n, k, i, "h", "", taille_zoom, spectral_lines, save, gamme, range, plot, mol_band)
            else:
                plot_lines(ax, f, l, path, synthetics, stardata, n, k, i, "k", "", taille_zoom, spectral_lines, save, gamme, range, plot, mol_band)
            # elif gamme == "visible":
            #     plot_lines(ax, f, l, path, synthetics, stardata, n, k, i, "visible", "", taille_zoom, spectral_lines, save, gamme, range, plot)
        
        save_dir = os.path.dirname(save)
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)
        plt.close()

def zoom_lines_analyse(ax, path, synthetics, molecules, stardata, j, m, k, i) : 
    # AX = syntspec(path+filename)

    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + m), stardata.get("v_"+ m)),stardata.get("flux_" + m),k)

    ax[0].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=1, color='black', label="Spectre observé : " + stardata.get("starname"))
    ax[1].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=1, color='black', label="Spectre observé : " + stardata.get("starname"))
    # ax[2].plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=0.8, color='lightgray', label="Spectre observé : " + stardata.get("starname"))

    for mol in molecules:
        AXb = syntspec(path + mol + j)
        ax[0].plot(AXb['wavelen'], AXb['flux'], label="Synth : atom +" + molecules.get(mol), linewidth=1)

    for synth in synthetics : 
        AX2 = syntspec(path+synth)
        ax[1].plot(AX2['wavelen'], AX2['flux'], label = "Synth : "+ synthetics.get(synth), linewidth = 1)

    # ax[0].plot(AX['wavelen'], AX['flux'], label="Synth : "+ filename, linewidth=1)
    # ax[2].plot(AX2['wavelen'], AX2['flux'], label = "Synth : atom + mols", color = "royalblue", linewidth = 0.8)

    # ax[2].set_xlim(k-5,k+5)
    ax[1].set_xlim(k-10,k+10)
    ax[0].set_xlim(k-10,k+10)

    for ax in ax : 
        ax.set_xlabel("Longueur d'onde (Å)", fontsize  = 12)
        ax.set_ylabel("Flux normalisé", fontsize  = 12)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major')
        ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor')
        ax.set_ylim(0.2, 1.3)
        ax.axvline(x=k, ymin=0., ymax=0.9, color='black',  linewidth = 0.5)  # Ajoute la barre verticale à la position 'zoom'
        ax.text(k, 1.2, s = i, color='black', fontsize=10, ha='center')  # Ajoute le nom au-dessus de la barre
        ax.legend(loc = "lower left", fontsize = 10)


def lines_analyse(lines, path, synthetics, molecules, stardata) :
    for i in list(lines.keys()):
        for k in lines.get(i):
            f = plt.figure(figsize=(15,8), dpi = 400)
            gs = f.add_gridspec(2, hspace=0.2)
            ax = gs.subplots(sharex=False, sharey=True)
            if k < 18500 : 
                zoom_lines_analyse(ax, path,synthetics, molecules, stardata, "", "h", k, i)

            elif k < 22500 : 
                zoom_lines_analyse(ax, path, synthetics,molecules, stardata, "2", "k", k, i)
            
            else :  
                zoom_lines_analyse(ax, path,synthetics, molecules, stardata, "3", "k", k, i)

            # plt.savefig("/Users/margauxvandererven/Documents/unif2023-2024/spectre_IR/output/"+stardata.get("starname")+"/"+i+"_"+str(k)+".png")


def zoom_lines_analyse_simple(ax, path,lines_mol, synthetics, stardata, m, k, i) : 
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + m), stardata.get("v_"+ m)),stardata.get("flux_" + m),k)

    ax.plot(normal['z_wavelen'], normal['flux_normalised'], linewidth=0.8, color='lightgray', label="Spectre observé : " + stardata.get("starname"))

    for synth in synthetics : 
        AX2 = syntspec(path+synth)
        ax.plot(AX2['wavelen'], AX2['flux'], label = "Synth : "+ synth, linewidth = 0.8)

    ax.set_xlim(k-10,k+10)
    
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major')
    ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor')
    ax.set_ylim(0.2, 1.3)
    ax.axvline(x=k, ymin=0., ymax=0.9, color='gray',  linewidth = 0.5)  # Ajoute la barre verticale à la position 'zoom'
    ax.text(k, 1.2, s = i, color='gray', fontsize=10, ha='center')  # Ajoute le nom au-dessus de la barre
    ax.legend(loc = "lower left", fontsize = 6)

    for d in list(lines_mol.keys()):
        t = lines_mol.get(d)
        ax.axvline(x=t, ymin=0., ymax=0.9, color='gray',  linewidth = 0.5)  # Ajoute la barre verticale à la position 'zoom'
        # ax.text(t, 1.2, s = d, color='gray', fontsize=10, ha='center')  # Ajoute le nom au-dessus de la barre


def lines_analyse_simple(lines, path, lines_mol, synthetics, stardata) :
    for i in list(lines.keys()):
        for k in lines.get(i):
            f = plt.figure(figsize=(20,8), dpi = 400)
            gs = f.add_gridspec(1, hspace=0.2)
            ax = gs.subplots(sharex=False, sharey=True)
            if k < 18500 : 
                zoom_lines_analyse_simple(ax, path, lines_mol,synthetics, stardata, "h", k, i)
            elif k < 22500 : 
                zoom_lines_analyse_simple(ax, path, lines_mol, synthetics, stardata, "k", k, i) 
            else :  
                zoom_lines_analyse_simple(ax, path, lines_mol,synthetics, stardata, "k", k, i)