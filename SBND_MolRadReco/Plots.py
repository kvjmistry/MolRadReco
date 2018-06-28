"""
This code will be used to make plots of the data from the MolRadReco Module.
Longitudinal Profile
Peak Maximum
Moliere Radius

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
    return fTmax(E, E_c, X_0) + (0.08 * Z + 9.6 ) * X_0 

# Values are from showers averaged over 500 events

# Variable Definitions
Z = 18      # Atomic Number of Argon
X_0 = 14.1  # Radiation length

# Energy
Energy =  np.arange(200, 3200, 200)
E_s    =  0.511 * np.sqrt(137 * 4 * np.pi)  # Scale energy
E_c    =  35.2                              # Critical Energy

# The longitudinal profie in truth at 95% containment
Lon_Truth  = np.array([108, 116, 128,132, 136, 140, 144,144, 148, 148, 148, 152, 152, 152, 156]) + 51
Lon_Theory = fLongitudinal(Energy, E_c, Z, X_0)
Lon_Diff   = Lon_Theory - Lon_Truth

# The maximum of the longitudinal profile
TMax_Truth        = np.array([18, 26, 30, 34, 42, 42, 42, 42, 54, 42, 42, 54, 54, 54, 54])                                    # Peak Max values unfitted
TMax_Truth_Fit    = np.array([18, 22.4, 27.5, 33.6, 37.0, 39.7, 39.8, 42.9, 46.5, 45.9, 48.6, 48.2, 50.2, 50.8, 52.6  ] )     # Fitted Peak max values
TMax_Theory       = fTmax(Energy, E_c, X_0)                                                                                   # Get the TMax Variable from Theory
TMax_Diff         = TMax_Theory - TMax_Truth_Fit # Residuals

# Moliere Radius
R_m_3D     = np.array([14.3, 13.7, 14.4, 13.9, 13.6, 13.7, 13.7, 13.7, 13.7, 13.7, 14.0, 13.7, 13.6, 13.5, 13.7])     # 3D Truth R_m
R_m_3D_2   = np.array([23.1, 22.2, 23, 22.7, 22.4, 22.4, 22.4, 22.8, 22.3, 22.5, 23.6, 22.6, 22.4, 22.3, 22.4 ])      # 3D Truth 2 * R_m
R_m_PDG    = 8.5  # PDG and Amaldi Definition of Moliere Radius
R_m_Fabjan = 11   # Fabjan and Amaldi Definition of Moliere Radius
# ----------------------- Plot the TMax Graph --------------------

plt.figure(0)
plt.plot(Energy , TMax_Theory , 'r-',  markersize=3, label = 'Theory')
plt.plot(Energy , TMax_Truth_Fit , 'go',  markersize=3, label = 'Truth')
plt.ylabel("TMax [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() 
plt.savefig("./Plots/TMax.png", dpi=500)

plt.figure(1)
plt.plot(Energy , TMax_Diff , 'bo',  markersize=3, label = 'TMax')
plt.ylabel("Diff [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() 
plt.savefig("./Plots/TMaxDiff.png", dpi=500)

# ----------------------- Plot the L @ 95% Graph --------------------
plt.figure(2)
plt.plot(Energy , Lon_Theory , 'r-',  markersize=3, label = 'Theory')
plt.plot(Energy , Lon_Truth , 'go',  markersize=3, label = 'Truth')
plt.ylabel("Conatinment at 95% [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() 
plt.savefig("./Plots/Lon_95%.png", dpi=500)

plt.figure(3)
plt.plot(Energy , Lon_Diff , 'bo',  markersize=3, label = 'Lon(95%)')
plt.ylabel("Diff [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(fontsize = "small") 
plt.savefig("./Plots/Lon95%_Diff.png", dpi=500)

# ----------------------- Section 2: Moliere Radius --------------------

plt.figure(4)
plt.plot(Energy , R_m_3D , 'ro',  markersize=3, label = 'Truth')
plt.axhline(y=R_m_Fabjan, color='b', linestyle='--', label = "Fabjan+ Amaldi")
plt.axhline(y=R_m_PDG, color='k', linestyle='--', label = "PDG + Amaldi")
plt.ylabel("Moliere Radius [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize = 13.5)
plt.savefig("./Plots/Moliere_Radius_SBND.png", dpi=500)

plt.figure(5)
plt.plot(Energy , 2 * R_m_3D , 'go',  markersize=3, label = 'Truth (2 * 90% Containment)')
plt.plot(Energy , R_m_3D_2 , 'ro',  markersize=3, label = 'Truth (95% Containment)')
plt.axhline(y=2 * R_m_Fabjan, color='b', linestyle='--', label = "Fabjan + Amaldi")
plt.axhline(y=2 * R_m_PDG, color='k', linestyle='--', label = "PDG + Amaldi")
plt.ylabel("2 * Moliere Radius [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize = 13.5)
plt.savefig("./Plots/Moliere_Radius_2_SBND.png", dpi=500)
