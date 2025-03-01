from imports import *

f = plt.figure(figsize=(20, 5), dpi = 200)
gs = f.add_gridspec(1, hspace=0.2)
(ax) = gs.subplots(sharex=False, sharey=False)

ax.plot(stardata.get("wavelen_visible"), stardata.get("flux_visible"), linewidth=0.5, color = 'black', label = starname)

ax.set(xlabel="Longueur d'onde (Å)", ylabel="Flux")
ax.xaxis.set_tick_params(direction = 'in', length = 5, which = 'major')
ax.xaxis.set_tick_params(direction = 'in', length = 2, which = 'minor')
cid = f.canvas.mpl_connect('button_press_event', on_click)

plt.show()