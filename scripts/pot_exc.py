import pprint
from stardata_BD22 import *

raie_Fe = {
 14745.387: [14745, 14746],
 14814.734: [14814.25, 14815.5],
 14826.408: [14825.9, 14827],
 14897.405: [14896.9, 14897.9],
 14988.778: [14988.3, 14989.4],
 15017.7: [15017.25, 15018.17],
 15077.287: [15076.78, 15077.8],
 15194.49: [15194, 15195.0],
 15343.788: [15343.33, 15344.33],
 15394.673: [15394.2, 15395.2],
 15395.718: [15395.2, 15396.35],
 15591.49: [15590.95, 15592.0],
 15648.51: [15647.99, 15649.0],
 15723.586: [15723, 15724.1],
 15741.918: [15741.43, 15742.46],
 15818.142: [15817.58, 15818.67],
 15821.712: [15821.21, 15822.2],
 15822.817: [15822.3, 15823.38],
 15911.302: [15910.86, 15911.87],
 15964.865: [15964.25, 15965.5],
 16040.654: [16040.2, 16041.18],
 16042.716: [16042.17, 16043.27],
 16125.899: [16125.43, 16126.53],
 16165.029: [16164.4, 16165.64],
 16180.9: [16180.44, 16181.4],
 16198.503: [16197.85, 16199],
 16284.769: [16284.22, 16285.5],
 16316.32: [16315.76, 16317],
 16318.691: [16318.1, 16319.3],
 16324.452: [16323.95, 16325.1],
 16436.621: [16436.15, 16437.15],
 16466.922: [16466.24, 16467.63],
 16486.667: [16486.06, 16487.3],
 16506.293: [16505.75, 16506.82],
 16517.223: [16516.66, 16517.74],
 16645.874: [16645.3, 16646.46],
 17204.297: [17203.74, 17204.91],
 17706.615: [17705.95, 17707.2],
 17721.086: [17720.51, 17721.86],
 17721.373: [17720.5, 17721.82],
 17771.123: [17770.5, 17771.8],
 17930.16: [17929.56, 17930.93],
 17932.6: [17931.95, 17933.25]
 }

raies_atom_H = list(raie_Fe.keys())
raies_atom_K = [
     19791.865, 19923.344,20281.084,21238.467,22257.107, 22260.18,22392.879,22419.977,22832.363
]

data=[]
with open("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/Linelists/Sophie_IGRINS/9000-15000_10042024.bsyn","r",encoding="utf8") as file :
        lines = file.readlines()
        for line in lines:
            parts=line.split()
            if len(parts)>3 and "'Fe" in parts:
                # print(parts)           
                wavelength = float(parts[0])
                excitation_potential = float(parts[1])
                loggf = float(parts[2])
                ew =  10**(loggf - (5040/4000)*excitation_potential)
                if ew >0.000000000000000001 and wavelength>14500:
                    data.append((wavelength, excitation_potential))
                else:
                     pass
            else:
                 pass
            
with open("/Users/margauxvandererven/Library/CloudStorage/OneDrive-UniversitéLibredeBruxelles/memoire/Linelists/Sophie_IGRINS/turbo_atoms.20180901_TS2020_transitions_mod_xx_ABO.txt","r",encoding="utf8") as file :
        lines = file.readlines()
        for line in lines:
            parts=line.split()
            if len(parts)>3 and "'FE" in parts:
                # print(parts)           
                wavelength = float(parts[0])
                excitation_potential = float(parts[1])
                loggf = float(parts[2])
                ew =  10**(loggf - (5040/4000)*excitation_potential)
                if ew >0.00000000000000001:
                    data.append((wavelength, excitation_potential))
                else:
                     pass
            else:
                 pass

# with open("../Linelists/Sophie_IGRINS/17000-25000_10042024.bsyn","r",encoding="utf8") as file :
#         lines = file.readlines()
#         for line in lines:
#             parts=line.split()
#             if len(parts)>3 and "'Fe" in parts:
#                 # print(parts)           
#                 wavelength = float(parts[0])
#                 excitation_potential = float(parts[1])
#                 loggf = float(parts[2])
#                 ew =  10**(loggf - (5040/4000)*excitation_potential)
#                 if ew >0.0000001:
#                     data.append((wavelength, excitation_potential))
#                 else:
#                      pass
#             else:
#                  pass

excitation_dict = {item[0]: item[1] for item in data}

raies_Fe = {wl: excitation_dict.get(wl, "Non trouvé") for wl in raies_atom_H+raies_atom_K}

grand_pot_ = {wl: excitation_dict.get(wl, "Non trouvé")for wl in lines_BD22.get("Fe I")}
grand_pot_2 = {key: value for key, value in grand_pot_.items() if value !="Non trouvé"}
grand_pot = {key: value for key, value in excitation_dict.items() if value > 10}

# print(grand_pot)
# pprint.pprint(results) 