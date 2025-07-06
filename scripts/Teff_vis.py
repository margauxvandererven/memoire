from imports import *

raie_Fe_vis=[7389.398,7418.667,7443.022,7461.263,7498.530,7540.430,7568.899,7586.018,7832.196,7937.139,8108.320,8698.706,8699.454,8710.404,8729.144,8747.425,8763.966]

# test_raie_Fe=

path_to_synth = "/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/BD-221742/synth_margaux/"
synth={"4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_Feabu_-20.conv":"$\epsilon_{Fe}$=0",
        "4000g1.0z-0.50m1.0t02a+0.20c+0.346n+0.00o+0.20r+0.00s+0.00.mod_4000-6000_vis_Feabu_7.2.conv":"log$\epsilon_{Fe}$= 7.2"}
range_dict={}
for i in test_raie_Fe:
    if 5000<i<6000:
        print("Processing line at wavelength:", i)
        zoom_lines({"":[i]}, path_to_synth, synth,stardata,2,lines_BD22_vis,gamme="visible")
#     print(range_dict)

# {4939.687:[4939.40,4939.92], 5434.524:[5434.26, 5434.90], 7418.667: [np.float64(7418.453), np.float64(7418.947)], 7443.022: [np.float64(7442.77), np.float64(7443.437)], 7461.263: [np.float64(7461.27), np.float64(7461.852)], 7498.53: [np.float64(7498.271), np.float64(7498.873)], 7540.43: [np.float64(7540.162), np.float64(7540.717)], 7568.899: [np.float64(7568.559), np.float64(7569.2)], 7586.018: [np.float64(7585.583), np.float64(7586.37)], 7832.196: [np.float64(7831.765), np.float64(7832.627)], 7937.139: [np.float64(7936.676), np.float64(7937.454)]}

abu_fe=np.arange(6.7,7.7,0.1)
# print(abu_fe)
