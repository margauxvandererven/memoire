# %matplotlib qt
# %load_ext autoreload
# %autoreload 2
# %reload_ext autoreload
from imports import *
repertory_memoire="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/"


# minimum = find_peaks_element("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/",
#                    "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_Feabu_7.2.conv",
#                    "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_sans12CH.conv")["minima"]
# for i in minimum:
#     if 4500 <i < 4700:
#         path_to_synth = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
#         synth={"4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_Feabu_7.2.conv":"avec 12CH",
#                     "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_sans12CH.conv": "sans 12CH"}
#         zoom_lines({"":[i]}, path_to_synth, synth,stardata,2,lines_BD22_vis,gamme="visible")

path_to_synth = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
Cabu=[7.88, 7.9, 7.8]
synth={}
# synth["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_sans12CH.conv"]= "sans 12CH"
synth["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4600-4800_vis_sans_C2.conv"]= "sans 12C12C"

for abu in Cabu:
    synth[f"4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4200-4400_vis_Cabu_{abu}.conv"]=f"log $\epsilon_{{O}}$={abu}"

# zoom_lines({"":[4700]}, path_to_synth, synth,stardata,10,lines_BD22_vis,gamme="visible")

fig = plot_lines_interactive(path_to_synth, synth, stardata, 4700, 10, lines_BD22_vis, gamme="visible")
fig.show()

# 4302.42,4306.06
# 4310.69 4312.45
# 4322.76 4324.62