"""
This code will be used to make plots of the data from the MolRadReco Module.
Longitudinal Profile
Peak Maximum

"""

import numpy as np
import matplotlib.pyplot as plt

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

### --- ###### --- ###### --- ###
# --- ### Start code  ### --- ###
### --- ###### --- ###### --- ###

# Functions
# Calculated the maximum of the longitudinal profile
def fTmax(E, E_c, X_0):
    return X_0 * (np.log(E/E_c) - 1) # -1 for el, -0.5 for gamma

# Returns the 95% cintainment estimation in theory
def fLongitudinal(E, E_c, Z, X_0):
    return fTmax(E, E_c, X_0) + (0.08 * Z + 9.6* X_0  ) - 32


# Values are from showers averaged over 500 events

# Variable Definitions
Z = 18      # Atomic Number of Argon
X_0 = 14.1  # Radiation length

# Energy
Energy =  np.arange(200, 3200, 200)
E_s    =  0.511 * np.sqrt(137 * 4 * np.pi)  # Scale energy
E_c    =  35.2                              # Critical Energy

# The longitudinal profie in truth at 95% containment
Lon_Truth = np.array([108, 116, 128,132, 136, 140, 144,144, 148, 148, 148, 152, 152, 152, 156]) 
Lon_Theory = fLongitudinal(Energy, E_c, Z, X_0)

# The maximum of the longitudinal profile
TMax_Truth    = np.array([18, 26, 30, 34, 42, 42, 42, 42, 54, 42, 42, 54, 54, 54, 54])
TMax_Theory = fTmax(Energy, E_c, X_0) # Get the TMax Variable from Theory

# ----------------------- Plot the TMax Graph --------------------


plt.figure(0)
plt.plot(Energy , TMax_Theory , 'r-',  markersize=3, label = 'Theory')
plt.plot(Energy , TMax_Truth , 'go',  markersize=3, label = 'Truth')
plt.ylabel("TMax [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() # put legend at the top right of the plot
plt.savefig("./Plots/TMax.png", dpi=500)

# ----------------------- Plot the L @ 95% Graph --------------------
plt.figure(1)
plt.plot(Energy , Lon_Theory , 'r-',  markersize=3, label = 'Theory')
plt.plot(Energy , Lon_Truth , 'go',  markersize=3, label = 'Truth')
plt.ylabel("Conatinment at 95% [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() # put legend at the top right of the plot
plt.savefig("./Plots/Lon_95%.png", dpi=500)