# %matplotlib qt
# %load_ext autoreload
# %autoreload 2
# %reload_ext autoreload
from imports import *
repertory_memoire="/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/"


minimum = find_peaks_element("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/",
                   "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_6000-8000_vis_sansCN.conv",
                   "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_6000-8000_vis_tout.conv")["minima"]


path_to_synth = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
Cabu=[7.88, 7.9, 7.8]
synth={}
# synth["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_sans12CH.conv"]= "sans 12CH"
synth["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_6000-8000_vis_sansCN.conv"]= "sans CN"
synth["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_6000-8000_vis_tout.conv"]= "avec CN"
synth_={}
synth_["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_sansCN.conv"]="sans CN"
synth_["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_tout.conv"]= "avec CN"

synth_2={}
synth_2["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_3550-4050_vis_sansCN.conv"]="sans CN"
synth_2["4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_3550-4050_vis_tout.conv"]= "avec CN"

list=list(np.arange(3600, 4050, 100))
# list_=np.arange(5500, 6100, 100)

CN_band=[[3860,3985],[4210,4240],[6320,6330]]

# for i in list:
#     if i < 6000:
#         synth_1=synth_
#     else:
#         synth_1=synth

#     zoom_lines({"":[i]}, path_to_synth, synth_1,stardata,50,lines_BD22_vis,gamme="visible")

zoom_lines({"":list}, path_to_synth, synth_2, stardata, 50,lines_BD22_vis,gamme="visible", save=f"../output/mol_vis/CN_{list[0]}-{list[-1]}", mol_band=CN_band)


# print(syntspec(path_to_synth + "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_Cabu_-20.conv")["header"])

# for abu in Cabu:
#     synth[f"4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4200-4400_vis_Cabu_{abu}.conv"]=f"log $\epsilon_{{O}}$={abu}"

# zoom_lines({"":[6320]}, path_to_synth, synth,stardata,20,lines_BD22_vis,gamme="visible")

# fig = plot_lines_interactive(path_to_synth, synth, stardata, 4700, 10, lines_BD22_vis, gamme="visible")
# fig.show()

# 4302.42,4306.06
# 4310.69 4312.45
# 4322.76 4324.62


#bande de CN:
#