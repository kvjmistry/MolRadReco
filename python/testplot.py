import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from lmfit import Minimizer, Parameters, report_fit, Model # lmfit for different style of fitting.

# These can be overwritten in the plot commands
#plt.close("all")  # Close all open figures

# Set the default plot variables
plt.rcParams['axes.linewidth'] = 2                 # Sets the boarder width 
plt.rcParams.update({'errorbar.capsize': 1.2})    # Gives a cap to the errorbars
#plt.rcParams["font.weight"] = "bold"              # sets all font weights to bold
plt.rcParams["font.size"] = 20                    # Sets all font sizes to bold
#plt.rcParams["axes.labelweight"] = "bold"         # Sets all label fontweights to bold 
plt.rcParams["lines.markersize"] = 4
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams['figure.figsize'] = 8, 8 # Sets the figure size

#choose the figure quality
figqual = 100

# Fit box style
props = dict(boxstyle='round', fc="white")

def fTmax(x, E_c, X_0):
    return X_0 * (np.log(x/E_c) - 1) # -1 for el, -0.5 for gamma

def fTmax_Residual(params, E, y):
    
    X_0 = params['X0']
    E_c = params['Ec']

    model = X_0 * (np.log(E/E_c) - 1)

    return y - model





# Variable Definitions
Z = 18      # Atomic Number of Argon
X_0 = 14.1  # Radiation length

# Energy
Energy =  np.arange(200, 3200, 200)
E_s    =  0.511 * np.sqrt(137 * 4 * np.pi)  # Scale energy
E_c  =  35.2                                # Critical Energy

TMax_Truth        = np.array([18, 26, 30, 34, 42, 42, 42, 42, 54, 42, 42, 54, 54, 54, 54])   

params = Parameters()
params.add('X0', value = X_0, min = 0,   max = 30, vary = True )
params.add('Ec', value = E_c, min = 10,  max = 60,  vary = True)


# do fit, here with leastsq model
minner = Minimizer(fTmax_Residual, params, fcn_args=(Energy, TMax_Truth))
result = minner.minimize()

# calculate final result
final = TMax_Truth + result.residual

result.params.pretty_print()

print result.params['X0'].value

# write error report
#report_fit(result)