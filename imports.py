# Matplotlib configuration
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Data manipulation
import numpy as np
from scipy.interpolate import griddata, interp1d
from scipy.optimize import minimize

# Custom modules
from stardata_BD22 import *
from zoom_raies import *
from minimisation_chi_2 import *
import pot_exc
from scipy.signal import find_peaks

# Plotting style
from tueplots import fonts
plt.rcParams.update(fonts.neurips2021())

# Other utilities
from pprint import pprint