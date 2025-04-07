import numpy as np

def ratio_C(ratio):
    p_2 = 1 / (ratio+1)
    p_1 = ratio / (ratio+1)
    print("12C : ", p_1, " & 13C : ", p_2)

def ratio_O(ratio):
    p_2 = 1 / (ratio+1)
    p_1 = ratio / (ratio+1)
    print("16O : ", p_1, " & 17O : ", p_2)

# ratio_O(100)
ratio=input("ratio C: ")
ratio_C(np.float64(ratio))
# 12meilleur