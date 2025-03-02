import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.optimize import minimize
from sklearn.linear_model import LinearRegression

from stardata_BD22 import *
from zoom_raies import *
from minimisation_chi_2 import *
# import pot_exc
# from Teff_clean import *
from scipy.signal import find_peaks
from utilities import *
from tueplots import fonts
plt.rcParams.update(fonts.neurips2021())

from pprint import pprint