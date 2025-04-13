from imports import *
from scipy import stats

def on_click(event):
    if event.xdata is not None and event.ydata is not None:
        coords = f"{event.xdata:.2f}"
        print(coords)
        # Copier les coordonnées dans le presse-papiers
        import pyperclip
        pyperclip.copy(coords)

def analyse_chi2(raies, ABU, variable,stardata,spectral_lines,iteration,minimisation=None, abu_to_plot=None,plot=None, name=None,save=None):
    """
    raies : dictionnaire des raies {raie:[start, end]}
    ABU : liste des abondances à tester
    chi_final : dictionnaire pour stocker les résultats {raie: [start, end, abu, min chi2]}
    variable : paramètre libre à minimiser
    abu_to_plot : liste des abondances à plotter
    name : nom de la raie
    """
    chi_final={}
    path_to_synth="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
    for wavelength in raies:
        try:
            print(f"Processing wavelength: {wavelength}")
            start=raies.get(str(wavelength))[0]
            end=raies.get(str(wavelength))[1]

            wavelength=np.float64(wavelength)
            start_range = wavelength - 2
            end_range = wavelength + 2
            
            if name=="OH":
                def format_number(num):
                    str_num = f"{num:.3f}"
                    if str_num.endswith('00'): 
                        return f"{num:.1f}"
                    elif str_num.endswith('0'): 
                        return f"{num:.2f}"
                    return str_num
            elif name=="12C16O":
                def format_number(num):
                    str_num = f"{num:.2f}"
                    if str_num.endswith('0'): 
                        return f"{num:.1f}"
                    return str_num
            else:
                return f"Aucun fichier de synthèse pour {name}"
        
            range = f"{format_number(start_range)}-{format_number(end_range)}"
            # range="14630-22900"
            synth={}
            model='4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod'
            if minimisation is not None: 
                for abu in ABU:
                    try:
                        # filename = f"{model}_{range}_{variable}abu_{abu:.2f}_{round}.conv"
                        filename = f"{model}_{range}_{variable}abu_{abu:.2f}_{iteration}.conv"
                        if os.path.exists(os.path.join(path_to_synth, filename)):
                            synth[filename] = f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
                        else:
                            # filename = f"{model}_{range}_{variable}abu_{str(abu)}_{round}.conv"
                            filename = f"{model}_{range}_{variable}abu_{str(abu)}_{iteration}.conv"
                            if os.path.exists(os.path.join(path_to_synth, filename)):
                                synth[filename] = f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
                            else:
                                print(f"Fichier non trouvé pour abu={abu}: {filename}")
                                continue
                    except Exception as e:
                        print(f"Erreur avec abu={abu}: {str(e)}")
                        continue
                    # print(synth)
                chi_squared_values = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=start,end=end)["chi_squared_values"]
                dof=chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=start,end=end)["dof"]
                chi_minimisation_ABU(ABU, chi_squared_values, variable,wavelength, name, chi_final,dof=dof, plot=plot, 
                                    #  save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/"+str(wavelength)+"/"+str(wavelength)+"_minimisation"
                                    #  save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final_Fe/"+str(wavelength)+"/"+str(wavelength)+"_minimisation"
                                     )

            if abu_to_plot is not None:
                synth_plot={}
                for abu2 in abu_to_plot:
                    try:
                        # filename = f"{model}_{range}_{variable}abu_{abu2:.2f}_{round}.conv"
                        filename = f"{model}_{range}_{variable}abu_{abu2:.2f}_{iteration}.conv"
                        if os.path.exists(os.path.join(path_to_synth, filename)):
                            synth_plot[filename] = f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                        else:
                            # filename = f"{model}_{range}_{variable}abu_{str(abu2)}_{round}.conv"
                            filename = f"{model}_{range}_{variable}abu_{str(abu2)}_{iteration}.conv"
                            if os.path.exists(os.path.join(path_to_synth, filename)):
                                synth_plot[filename] = f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                            else:
                                print(f"Fichier non trouvé pour abu={abu2}: {filename}")
                                continue
                    except Exception as e:
                        print(f"Erreur avec abu={abu2}: {str(e)}")
                        continue
                # synth_plot["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_"+range+"_sans_"+name+".conv"]= "sans "+name
                if wavelength < 18500:
                    synth_plot["../../syntspec/BD-221742b/Fe/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_14500-18500_BD-221742_Feabu_-20.conv"]= "sans Fe"
                else:
                    synth_plot["../../syntspec/BD-221742b/Fe/4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_19750-22900_BD-221742_Feabu_-20.conv"]= "sans Fe"
                plot_zone_chi2_simple(wavelength, path_to_synth, synth_plot, stardata, spectral_lines, size_police=14,size_trace=(1., 8),name=name, start=start,end=end,
                                    # save="/Users/margauxvandererven/OneDrive  - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/"+str(wavelength)+"/"+str(wavelength)+"_zone", plot=plot
                                    # save="/Users/margauxvandererven/OneDrive  - Université Libre de Bruxelles/memoire/output/final_Fe/"+str(wavelength)+"/"+str(wavelength)+"_zone", plot=plot
                                    )
        except Exception as e: 
            print(f"Error processing wavelength {wavelength}: {str(e)}")
            continue
    # print(chi_final)
    # if minimisation is not None:
    #     abu_plot(chi_final,variable,size_police=20,
    #             #  save="/Users/margauxvandererven/OneDrive - Université Libre de Bruxelles/memoire/output/final/"+name+"/"+round+"/Oabu_"+round,
    #             #  save="../rédaction/images/plot_abu/"+name+"_final_"+round+".pdf"
    #              )
   
    # chi_final[k] = [start, end, min_log_e, error_minus, error_plus, min_chi_squared, r_squared, p_value]
    keys_to_remove = []
    for k, v in chi_final.items():
        if len(v) != 8:
            keys_to_remove.append(k)
    
    for k in keys_to_remove:
        del chi_final[k]

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

def analyse_chi2_CO(raies, ABU, variable,round,stardata,spectral_lines,minimisation=None, abu_to_plot=None, name=None,save=None):
    """
    raies : dictionnaire des raies {raie:[start, end]}
    ABU : liste des abondances à tester
    chi_final : dictionnaire pour stocker les résultats {raie: [start, end, abu, min chi2]}
    variable : paramètre libre à minimiser
    abu_to_plot : liste des abondances à plotter
    name : nom de la raie
    """
    chi_final={}
    path_to_synth="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
    for wavelength in raies:
        start=raies.get(str(wavelength))[0]
        end=raies.get(str(wavelength))[1]

        wavelength=np.float64(wavelength)
        if wavelength < 18500:
            range="15525-17275"
        else: 
            range="23050-24800"
        synth={}
        for abu in ABU:
            synth["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_"+range+"_"+variable+"abu_"+str(abu)+"_"+round+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu)}"
        if abu_to_plot:
            synth_plot={}
            for abu2 in abu_to_plot:
                synth_plot["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_"+range+"_"+variable+"abu_"+str(abu2)+"_"+round+".conv"]= f"log$\\epsilon_{{{variable}}}$ = {str(abu2)}"
                if wavelength < 18500:
                    synth_plot["../../syntspec/BD-221742b/sans_CO_mol"]= "sans "+name
                else :
                    synth_plot["../../syntspec/BD-221742b/sans_CO_mol2"]= "sans "+name

            plot_zone_chi2(wavelength, path_to_synth, synth_plot, stardata, spectral_lines,axes=(True,True), size_police=20,size_trace=(1.8, 10),name=name, start=start,end=end)
        
        if minimisation is not None:
            chi_squared_values = chi_2(path_to_synth, synth, stardata, wavelength, spectral_lines,chi_final,name=name,start=start,end=end)["chi_squared_values"]
            chi_minimisation_ABU(ABU, chi_squared_values, variable,wavelength, name, chi_final)
    
    if minimisation is not None:
        abu_plot(chi_final,variable,size_police=20,save=variable+"_abu_"+round)
        if save:
            with open(save, "w") as fichier:
                json.dump(chi_final, fichier, indent=4, ensure_ascii=False)

def plot_zone_chi2(k, path, synthetics, stardata, spectral_lines, axes=None, size_police=None,size_trace=(None,None), save=None, start=None, end=None, name=None):
    if k < 18500:
        j="h"
    else :
        j="k"
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + j), stardata.get("v_" + j)), stardata.get("flux_" + j), k)
    
    f = plt.figure(figsize=(12,10))
    gs = f.add_gridspec(2, hspace=0.2)
    ax = gs.subplots(sharex=False, sharey=True)
    cid = f.canvas.mpl_connect('button_press_event', on_click)
    colors = ['#3e92f7',  '#faaa2e','#db2f2f','#a39cf3', '#92c124']
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
    line_styles = [
    {'linestyle': '-', 'linewidth': 2},      # ligne continue épaisse
    {'linestyle': '--', 'linewidth': 2},     # tirets épais
    {'linestyle': '-.', 'linewidth': 2},     # point-tiret épais
    {'linestyle': ':', 'linewidth': 2.5},    # pointillés épais
    {'linestyle': '-', 'linewidth': 1.5}     # ligne continue fine
    ]
    ax[0].set_xlim(k-10,k+10)
    ax[1].set_xlim(k-2,k+2)
    ax[0].set_ylim(0.0, 1.4)
    ax[1].set_ylim(0.2, 1.4)
    if size_police is None:
        size_police = 10
    if size_trace[0] is not None:
        linewidth_trace=size_trace[0]
    else:
        linewidth_trace=2.4
    if size_trace[1] is not None:
        marker_size=size_trace[1]
    else:
        marker_size=14

    
    if axes == (True, None):
        ax[1].set_xlabel("$\lambda$ [Å]", fontsize  = size_police)
    elif axes == (None, True):
        ax[0].set_ylabel("F$_{\mathrm{norm}}$", fontsize  = size_police)
        ax[1].set_ylabel("F$_{\mathrm{norm}}$", fontsize  = size_police)
    elif axes == (True, True):
        ax[0].set_ylabel("F$_{\mathrm{norm}}$", fontsize  = size_police)
        ax[1].set_xlabel("$\lambda$ [Å]", fontsize  = size_police)
        ax[1].set_ylabel("F$_{\mathrm{norm}}$", fontsize  = size_police)
    else:
        pass

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

    raies_filtrees=[]
    for i, (z, element) in enumerate(all_wavelengths):
        closest_index = (np.abs(np.array(normal['z_wavelen']) - z)).argmin()
        if normal['flux_normalised'][closest_index] <= 0.9:
            raies_filtrees.append((z, element))
        else:
            pass

    for i, (z, element) in enumerate(raies_filtrees):
        if i > 0 and abs(z - raies_filtrees[i - 1][0]) < 0.7:
            text_height += 0.07  
        else:
            text_height = text_height_base  

        if k - 10 <= z <= k + 10:
            ax[0].axvline(x=z, ymin=0.75, ymax=0.8, color='black', linewidth=1.)
            ax[0].text(z, text_height, s=element, color='black', fontsize=size_police-5, ha='center')
            # print(normal['flux_normalised'][closest_index])
        if k - 2 <= z <= k + 2:
            ax[1].axvline(x=z,   ymin=0.75, ymax=0.8, color='black', linewidth=1.)
            ax[1].text(z, text_height, s=element, color='black', fontsize=size_police-5, ha='center')
    
    ax[0].scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=marker_size, 
        facecolors='none', 
        color='black', 
        linewidths=1.5
    #    label="Spectre observé : " + stardata.get("starname")
        )
    ax[1].scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=marker_size+4, 
        facecolors='none', 
        color='black',
        linewidths=2 
    #    label="Spectre observé : " + stardata.get("starname")
        )
    
    for i, synth in enumerate(synthetics) : 
        AX2 = syntspec(path+synth)
        ax[0].plot(AX2['wavelen'], AX2['flux'], label=synthetics.get(synth), color=colors[i % len(colors)], linewidth=linewidth_trace
        # **line_styles[i % len(line_styles)],
        )
        ax[1].plot(AX2['wavelen'], AX2['flux'], label=synthetics.get(synth), color=colors[i % len(colors)], linewidth=linewidth_trace+0.6
        # **line_styles[i % len(line_styles)],
        )

    for ax_ in ax :
        ax_.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
        ax_.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
        ax_.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
        ax_.xaxis.get_major_formatter().set_scientific(False)
        ax_.xaxis.get_major_formatter().set_useOffset(False)
        ax_.tick_params(axis = 'both', labelsize = size_police)
        if start is not None and end is not None:
            ax_.axvspan(start, end, color='lightgray', alpha=0.15)
            ax_.axvline(x=start, ymin=0, ymax=1.4, color='gray', linewidth=1)
            ax_.axvline(x=end, ymin=0, ymax=1.4, color='gray', linewidth=1)
        ax_.tick_params(axis='x', pad=11)
    ax[0].legend(fontsize = size_police, framealpha=0.8, facecolor='white', markerscale=0.2,edgecolor='white',numpoints=5, loc='lower left')
    # plt.tight_layout()
    if save is not None:
        plt.savefig(save + '.pdf', dpi=400, transparent=True, bbox_inches
                    ='tight')
    plt.show()
        # plt.savefig(save + '.pgf', backend='pgf')


def plot_zone_chi2_simple(k, path, synthetics, stardata, spectral_lines, size_police=None,size_trace=(None,None), save=None, start=None, end=None, name=None, plot=None):
    if k < 18500:
        j="h"
    else :
        j="k"
    normal = normalisation(redshift_wavelen(stardata.get("wavelen_" + j), stardata.get("v_" + j)), stardata.get("flux_" + j), k)
    
    f = plt.figure(figsize=(8,4))
    gs = f.add_gridspec(1, hspace=0.2)
    ax = gs.subplots(sharex=False, sharey=True)
    cid = f.canvas.mpl_connect('button_press_event', on_click)
    colors = ['#3e92f7',  '#faaa2e','#db2f2f','#a39cf3', '#92c124']
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
    line_styles = [
    {'linestyle': '-', 'linewidth': 2},      # ligne continue épaisse
    {'linestyle': '--', 'linewidth': 2},     # tirets épais
    {'linestyle': '-.', 'linewidth': 2},     # point-tiret épais
    {'linestyle': ':', 'linewidth': 2.5},    # pointillés épais
    {'linestyle': '-', 'linewidth': 1.5}     # ligne continue fine
    ]
    ax.set_xlim(k-1.5,k+1.5)
    ax.set_ylim(0.2, 1.4)
    if size_police is None:
        size_police = 10
    if size_trace[0] is not None:
        linewidth_trace=size_trace[0]
    else:
        linewidth_trace=2.4
    if size_trace[1] is not None:
        marker_size=size_trace[1]
    else:
        marker_size=14

    ax.set_xlabel("$\lambda$ [Å]", fontsize  = size_police)
    ax.set_ylabel("F$_{\mathrm{norm}}$", fontsize  = size_police)

    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    all_wavelengths = []
    for element, wavelengths in spectral_lines.items():
        for wavelength in wavelengths:
            all_wavelengths.append((wavelength, element))

    all_wavelengths.sort()  
  
    text_height_base = 1.2
    text_height = text_height_base  

    raies_filtrees=[]
    for i, (z, element) in enumerate(all_wavelengths):
        closest_index = (np.abs(np.array(normal['z_wavelen']) - z)).argmin()
        if normal['flux_normalised'][closest_index] <= 0.9:
            raies_filtrees.append((z, element))
        else:
            pass

    for i, (z, element) in enumerate(raies_filtrees):
        if i > 0 and abs(z - raies_filtrees[i - 1][0]) < 0.7:
            text_height += 0.07  
        else:
            text_height = text_height_base  

        if k - 2 <= z <= k + 2:
            ax.axvline(x=z,   ymin=0.75, ymax=0.8, color='black', linewidth=1.)
            ax.text(z, text_height, s=element, color='black', fontsize=size_police-5, ha='center')
    
    ax.scatter(normal['z_wavelen'], normal['flux_normalised'], marker='o',s=marker_size+4, 
        facecolors='none', 
        color='black',
        linewidths=2 
    #    label="Spectre observé : " + stardata.get("starname")
        )
    
    for i, synth in enumerate(synthetics) : 
        AX2 = syntspec(path+synth)
        ax.plot(AX2['wavelen'], AX2['flux'], label=synthetics.get(synth), color=colors[i % len(colors)], linewidth=linewidth_trace+0.6
        # **line_styles[i % len(line_styles)],
        )
    ax.xaxis.set_tick_params(direction = 'in', length = 10, which = 'major', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction = 'in', length = 10, which = 'major',top=True, bottom=True)
    ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_major_formatter().set_useOffset(False)
    ax.tick_params(axis = 'both', labelsize = size_police)
    if start is not None and end is not None:
        ax.axvspan(start, end, color='lightgray', alpha=0.15)
        ax.axvline(x=start, ymin=0, ymax=1.4, color='gray', linewidth=1)
        ax.axvline(x=end, ymin=0, ymax=1.4, color='gray', linewidth=1)
    ax.tick_params(axis='x', pad=11)
    ax.legend(fontsize = size_police, framealpha=0.8, facecolor='white', markerscale=0.2,edgecolor='white',numpoints=5, loc='lower left')
    # plt.tight_layout()
    # if save is not None:
    save_dir = os.path.dirname(save)
    if save_dir and not os.path.exists(save_dir):
        os.makedirs(save_dir)
    plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)

    if plot:
        plt.show()
    else:
        plt.close()


def get_nearest(x_query, x_vals):
    nearest_index = np.argmin(np.abs(x_vals - x_query))
    return nearest_index


def chi_squared(synth, observed, errors=None):
    errors=0.01
    return np.sum(((observed - synth)** 2 /errors**2 ))


def chi_2(path, synthetics, stardata, k, spectral_lines,chi_final, start, end, name=None, plot=None, save=None, size_police=None, axes=(None,None)):
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

    # chi_final[k]=[start, end]
    chi_final[k]=[wavelength_observed_filtered[0],wavelength_observed_filtered[-1]]

    synthetic_spectra = []
    observed_spectra_interpolated =[]
    observed_spectra_w_interpolated =[]
    synthetic_spectra_w = []
    chi_squared_values = []

    for syntha in synthetics:
            # print(syntha)
            synt_flux = np.array(zoom_syntspec(path, syntha, k, 10)["synt_flux"])
            synt_w = np.array(zoom_syntspec(path, syntha, k, 10)["synt_wavelen"])

            mask_synthetic = (synt_w >= start) & (synt_w <= end)

            wavelength_synthetic_filtered = synt_w[mask_synthetic]
            flux_synthetic_filtered = synt_flux[mask_synthetic]

            interpolator = interp1d(synt_w, synt_flux, kind='linear', fill_value="extrapolate")
            flux_synthetic_interpolated = interpolator(wavelength_observed_filtered)

            chi2 = chi_squared(flux_synthetic_interpolated, flux_observed_filtered)
            # chi2_reduce = chi2 / (len(flux_observed_filtered)-1)
            chi_squared_values.append(chi2)

            synthetic_spectra.append(flux_synthetic_interpolated)
            synthetic_spectra_w.append(wavelength_observed_filtered)

    if plot is True:
        plot_zone_chi2(k, path, synthetics, normal, spectral_lines, size_police=size_police,
                        save=save, start=wavelength_observed_filtered[0], end=wavelength_observed_filtered[-1], name=name, axes=axes)

    return {"chi_squared_values":chi_squared_values, "dof":len(flux_observed_filtered)-1}


def quadratic(x, a, b, c):
    return a * x**2 + b * x + c

def chi_minimisation_ABU(abu, chi_squared, element, k, raie, chi_final, dof, plot=None, save=None):
    start=chi_final[k][0]
    end=chi_final[k][1]

    abu = np.array(abu)
    chi_squared_red = np.array(chi_squared)/dof
    
    def chi2_func(x):
        interp = interp1d(abu, chi_squared_red, kind='quadratic')
        return interp(x)
    
    try:
        result = minimize(chi2_func, 
                        x0=abu[np.argmin(chi_squared_red)],
                        bounds=[(abu.min(), abu.max())],
                        method='Nelder-Mead')
        
        min_log_e = result.x[0]
        min_chi_squared = result.fun
        
        x_fine = np.linspace(abu.min(), abu.max(), 1000)
        interp = interp1d(abu, chi_squared_red, kind='quadratic')
        y_fine = interp(x_fine)

        y_pred = interp(abu)
        residuals = chi_squared_red - y_pred
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((chi_squared_red - np.mean(chi_squared_red))**2)
        r_squared = 1 - (ss_res / ss_tot)
        p_value = 1 - stats.chi2.cdf(min_chi_squared * dof, dof)
        
        mask = y_fine <= min_chi_squared + 1/dof
        x_masked = x_fine[mask]
        error_minus = min_log_e - x_masked[0]
        error_plus = x_masked[-1] - min_log_e
        
        
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(abu, chi_squared_red, color='darkblue', marker='x', 
                    label=f"Raie de {raie} en {k} Å")
        
        ax.plot(x_fine, y_fine, 'gray', alpha=0.5)
        
        ax.hlines(min_chi_squared+1/dof,x_masked[0]-0.01,x_masked[-1]+0.01,color='darkseagreen', linestyle='-', alpha=0.8,
                    label=r'$\chi^2_{min} + 1/dof$')
        ax.axvline(x=x_masked[0], color='darkseagreen', linestyle='--', alpha=0.8)
        ax.axvline(x=x_masked[-1], color='darkseagreen', linestyle='--', alpha=0.8)
        
        ax.scatter(min_log_e, min_chi_squared, color='red', marker='x',
                    label=f"Minimum: {min_log_e:.3f}$\pm^{{{error_plus:.3f}}}_{{{error_minus:.3f}}}$")

        # plt.legend()
        ax.xaxis.set_tick_params(direction = 'in', length = 5, which = 'major', top=True, bottom=True)
        ax.yaxis.set_tick_params(direction = 'in', length = 5, which = 'major',top=True, bottom=True)
        ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
        ax.tick_params(axis = 'both', labelsize = 16)
        ax.set_xlabel(f"$\\log \\epsilon_{{\\mathrm{{{element}}}}}$", fontsize=16)
        ax.set_ylabel(r"$\chi^2_{red}$", fontsize=16)
        ax.text(0.4, 0.95, f'$\\log \\epsilon_{{\\mathrm{{{element}}}}}$ = {min_log_e:.3f}$\pm^{{{error_plus:.3f}}}_{{{error_minus:.3f}}}$', 
            transform=ax.transAxes, fontsize=14,
            verticalalignment='top')
        ax.text(0.4, 0.9, f'$\chi^2_{{red,min}}$ = {min_chi_squared:.3f}', 
            transform=ax.transAxes, fontsize=14,
            verticalalignment='top')
        
        if save:
            save_dir = os.path.dirname(save)
            if save_dir and not os.path.exists(save_dir):
                os.makedirs(save_dir)
            plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)
        
        if plot:
            plt.show()
        else:
            plt.close()
            
        chi_final[k] = [start, end, min_log_e, error_minus, error_plus, min_chi_squared, r_squared, p_value]
        
    except Exception as e:
        print(f"Erreur lors de la minimisation pour la raie {k}: {str(e)}")

def chi_minimisation_ABU_ancien(abu, chi_squared, element, k, raie, chi_final, dof, plot=None, save=None):
    abu = np.array(abu)
    chi_squared_red = np.array(chi_squared)/dof
    params, _ = curve_fit(quadratic, abu, chi_squared_red)
    a, b, c = params  
    fit_x = np.linspace(abu.min(), abu.max(), 100)
    fit_y = quadratic(fit_x, a, b, c)

    min_log_e = -b / (2 * a)  
    min_chi_squared = quadratic(min_log_e, a, b, c)

    residuals = chi_squared_red - quadratic(abu, a, b, c)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((chi_squared_red - np.mean(chi_squared_red))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f"R² de l'ajustement : {r_squared:.3f}")
    
    p_value = 1 - stats.chi2.cdf(min_chi_squared*dof, dof)
    print(f"p-value : {p_value:.3f}")
    if p_value < 0.05:
        print("p-value < 0.05 : ajustement potentiellement problématique")
        if min_chi_squared > 3:
            print("Chi² réduit > 3 : mauvais ajustement")
        elif min_chi_squared < 0.1:
            print("Chi² réduit < 0.1 : possible sur-ajustement")

    if plot is not None:
        f = plt.figure(figsize=(8, 6))
        gs = f.add_gridspec(1, hspace=0.2)
        ax = gs.subplots(sharex=False, sharey=True)
        ax.scatter(abu, chi_squared_red, color="darkblue", marker="x", label=f"Raie de {raie} en {k} Å")  # Data points
        ax.plot(fit_x, fit_y, color="lightgray")  # Quadratic fit line
        ax.scatter(min_log_e, min_chi_squared, color="red", marker="x", label=f"Minimum en {min_log_e:.2f}")  # Minimum point
        ax.xaxis.set_tick_params(direction = 'in', length = 5, which = 'major', top=True, bottom=True)
        ax.yaxis.set_tick_params(direction = 'in', length = 5, which = 'major',top=True, bottom=True)
        ax.xaxis.set_tick_params(direction = 'in', length = 6, which = 'minor',top=True, bottom=True)
        ax.tick_params(axis = 'both', labelsize = 16)
        # Labeling
        ax.set_xlabel(f"$\\log \\epsilon_{{\\mathrm{{{element}}}}}$", fontsize=16)
        # ax.set_xlabel(element, fontsize=16)
        ax.set_ylabel(r"$\chi^2_{red}$", fontsize=16)
        discriminant = b**2 - 4*a*(c - (min_chi_squared + 1/dof))
        if discriminant >= 0:
            x1 = (-b - np.sqrt(discriminant))/(2*a)
            x2 = (-b + np.sqrt(discriminant))/(2*a)
            
            if plot is not None:
                # ...existing plotting code...
                # Ajouter les lignes verticales aux intersections
                ax.axvline(x=x1, color='darkseagreen', linestyle='--', alpha=0.5)
                ax.axvline(x=x2, color='darkseagreen', linestyle='--', alpha=0.5)
                ax.hlines(min_chi_squared+1/dof, x1-0.01, x2+0.01, color="darkseagreen", label="Minimum + 1/dof")

                # Ajouter le texte pour l'incertitude
                # ax.text(0.05, 0.95, f'Δ = {abs(x2-x1):.2f}', 
                #     transform=ax.transAxes, fontsize=14,
                #     verticalalignment='top')
                half_width = abs(x2-x1)/2
                # ax.text(0.05, 0.95, f'σ = {half_width}', 
                #     transform=ax.transAxes, fontsize=14,
                #     verticalalignment='top')
                ax.text(0.4, 0.95, f'log $\epsilon_{element}$ = {min_log_e:.3f} $\pm$ {half_width:.3f}', 
                    transform=ax.transAxes, fontsize=14,
                    verticalalignment='top')
                # ax.text(0.50, min_chi_squared+1/dof+0.01, f'$\chi^2_{{red,min}}$ + 1/dof', 
                #     transform=ax.transAxes, fontsize=14,color="darkseagreen",
                #     verticalalignment='top')
                ax.text(0.4, 0.9, f'$\chi^2_{{red,min}}$ = {min_chi_squared:.3f}', 
                    transform=ax.transAxes, fontsize=14,
                    verticalalignment='top')
        # plt.legend(loc="upper left", fontsize=16)
    if save is not None:
        # Extraire le chemin du dossier depuis le chemin complet
        save_dir = os.path.dirname(save)
        
        # Créer le dossier s'il n'existe pas
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
        plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)
    # plt.close()
    plt.show()
    # plt.show(block=False)
    # plt.pause(0.1)

    chi_final.get(k).append(min_log_e)
    chi_final.get(k).append(half_width)
    chi_final.get(k).append(min_chi_squared)
    chi_final.get(k).append(r_squared)
    chi_final.get(k).append(p_value)

def interpolated_chi_squared(x, data1, data2, chi_squared_values):
    ab, mt = x
    return griddata((data1, data2), chi_squared_values, (ab, mt), method='cubic')

# Fonction de minimisation et visualisation
def double_chi(abundances, macroturbulences, chi_squared_values,line, data, size_police=None, save=None):
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

    data[line]=[best_abundance
    , best_macroturbulence, min_chi_squared_value]

    print("Meilleure abondance estimée :", best_abundance)
    print("Meilleur ratio C/O estimée :", best_macroturbulence)
    print("Valeur minimale de χ² estimée :", min_chi_squared_value)

    if size_police is None:
        size_police = 10
    # Visualisation
    plt.figure(figsize=(10, 8))
    plt.contourf(grid_ab, grid_mt, grid_chi_squared, levels=250, cmap='RdYlBu_r', alpha=0.8,                # Étendre les couleurs
    antialiased=True)
    cbar=plt.colorbar(label='$\chi^2$')
    plt.scatter(abundances, macroturbulences, c='white',marker='x')
    plt.plot(best_abundance, best_macroturbulence, color='red',marker='+', markersize=12)
    plt.xlabel("log $\\epsilon_{\\mathrm{O}}$", fontsize=size_police)
    plt.ylabel("C/O", fontsize=size_police)
    plt.tick_params(labelsize = size_police)
    cbar.ax.tick_params(labelsize=size_police)
    cbar.set_label('$\chi^2$', fontsize=size_police)
    # plt.legend()
    plt.show()
    if save is not None:
        plt.savefig("rédaction/images/chi2/"+save+".pdf", dpi=400, transparent=True, bbox_inches='tight')


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

def abu_plot(raies_element, element_abu, save=None, size_police=None):
    x_vals=[]
    y_vals=[]
    errors_minus = [] 
    errors_plus = [] 
    for key in raies_element.keys():
        val=raies_element[key]
        p_value=np.float64(val[-1])
        if 0.05<p_value<0.95 and val[3]>0 and val[4]>0:
            x_vals.append(np.float64(key))
            y_vals.append(val[2])
            errors_minus.append(val[3])  
            errors_plus.append(val[4]) 

    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)
    errors = np.array([errors_minus, errors_plus])

    f = plt.figure(figsize=(10, 7))
    gs = f.add_gridspec(1)
    ax = gs.subplots(sharex=False, sharey=True)
    if size_police is None:
        size_police = 10
    # Scatter plot for raies
    ax.errorbar(x_vals, y_vals, yerr=errors, color="darkblue",capsize=5,  fmt='o')
    ax.set_xlabel("$\\lambda$ [Å]", fontsize=size_police)
    ax.set_ylabel(f"log $\\epsilon_{{\\mathrm{{{element_abu}}}}}$", fontsize=size_police)
    
    # Set tick parameters
    ax.xaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    ax.xaxis.set_tick_params(direction='in', length=5, which='minor', top=True, bottom=True)
    ax.yaxis.set_tick_params(direction='in', length=8, which='major', top=True, bottom=True)
    ax.tick_params(axis = 'both', labelsize = size_police)
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
    # print("Mean of filtered points:", filtered_mean)
    # print("Standard deviation of filtered points:", filtered_std)
    
    # Plot horizontal line for the filtered mean
    ax.hlines(y=y_mean, xmin=14500, xmax=23500, color="indianred", linestyle="--", linewidth=1.8, label=f"$\mu$ = {y_mean:.2f}")

    # Plot original shaded error band
    ax.fill_between(
        [14500, 23500], 
        y_mean - y_std, 
        y_mean + y_std, 
        color="gray", 
        alpha=0.1, 
        label=f"$\sigma$ = ±{y_std:.2f}"
    )
    ax.set_xlim(14500, 23500)
    # Display legend
    # plt.legend(fontsize=size_police)
    plt.text(0.4, 0.95, f"log $\\epsilon_{{\\mathrm{{{element_abu}}}}}$ = {y_mean:.2f} $\pm$ {y_std:.2f}", 
    transform=ax.transAxes, fontsize=size_police,
    verticalalignment='top')

    if save is not None:
        plt.savefig(save + ".pdf", dpi=600, bbox_inches='tight', transparent=True)
        
    plt.show()
    # plt.close()
    # Return filtered mean and standard deviation
    # return filtered_mean, filtered_std