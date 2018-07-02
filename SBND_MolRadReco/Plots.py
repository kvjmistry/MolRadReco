"""
This code will be used to make plots of the data from the MolRadReco Module.
Longitudinal Profile
Peak Maximum
Moliere Radius

"""

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
plt.rcParams["font.size"] = 18                    # Sets all font sizes to bold
#plt.rcParams["axes.labelweight"] = "bold"         # Sets all label fontweights to bold 
plt.rcParams["lines.markersize"] = 4
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams['figure.figsize'] = 8, 8 # Sets the figure size

# choose the figure quality
figqual = 100

# Fit box style
props = dict(boxstyle='square', fc="white", alpha = 1)

### --- ###### --- ###### --- ###
# --- ### Start code  ### --- ###
### --- ###### --- ###### --- ###

# Functions
# Calculated the maximum of the longitudinal profile
def fTmax(E, E_c, X_0):
    return X_0 * (np.log(E/E_c) - 1) # -1 for el, -0.5 for gamma

def fTmax_Residual(params, E, y):   
    X_0 = params['X_0']
    E_c = params['E_c']

    model = X_0 * (np.log(E/E_c) - 1)

    return (y - model)

# Returns the 95% cintainment estimation in theory
def fLongitudinal(E, E_c, Z, X_0):
    return fTmax(E, E_c, X_0) + (0.08 * Z + 9.6 ) * X_0 

def fLongitudinal_Residual(params, E, y):
    X_0 = params['X_0']
    E_c = params['E_c']
    Z   = params['Z']

    model = (np.log(E/E_c) - 1   + 0.08 * Z  + 9.6 ) * X_0

    return (y - model)

# Values are from showers averaged over 500 events

# Variable Definitions
Z = 18      # Atomic Number of Argon
X_0 = 14.1  # Radiation length

# Energy
Energy =  np.arange(200, 3200, 200)
E_s    =  0.511 * np.sqrt(137 * 4 * np.pi)  # Scale energy
E_c    =  35.2                                # Critical Energy

# =+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+
# =+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+

# The longitudinal profie in truth at 95% containment
Lon_Truth  = np.array([108, 116, 128,132, 136, 140, 144,144, 148, 148, 148, 152, 152, 152, 156]) + 0 # Shift was 51
Lon_Theory = fLongitudinal(Energy, E_c, Z, X_0)
Lon_Ratio  = Lon_Theory / Lon_Truth

# create a set of Parameters
p_Lon = Parameters()
p_Lon.add('X_0', value = X_0, min = 0,   max = 20,  vary = True)
p_Lon.add('E_c', value = E_c, min = 20,  max = 40,  vary = True)
p_Lon.add('Z',   value = Z,  vary = False)


# Do fit, here with leastsq model
minner_Lon = Minimizer(fLongitudinal_Residual, p_Lon, fcn_args=(Energy, Lon_Truth))
fitresult_Lon = minner_Lon.minimize()

# calculate the fit residuals
residual_Lon = fitresult_Lon.residual

# write error report
print "\n\n=====  Lon (95%) FIT REPORT ===== "
report_fit(fitresult_Lon)
print "===== ===================== ===== \n\n"
 
X_0_fit = fitresult_Lon.params['X_0'].value
E_c_fit = fitresult_Lon.params['E_c'].value 

# Generate the fit to plot 
Lon_Truth_fit = fLongitudinal(Energy, E_c, Z, X_0_fit)

# Create a string for displaing on the plot
textstr_Lon = 'Fit Parameters \n---------------------\n $X_0=%.1f$ cm \n $ E_c=%.1f$ MeV \n $\chi_r^2=%.1f$ \n' % (X_0_fit, E_c_fit,fitresult_Lon.redchi )


# =+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+
# =+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+

# The maximum of the longitudinal profile (TMAX)
TMax_Truth        = np.array([18, 26, 30, 34, 42, 42, 42, 42, 54, 42, 42, 54, 54, 54, 54])                                    # Peak Max values unfitted
TMax_Truth_PFit   = np.array([18, 22.4, 27.5, 33.6, 37.0, 39.7, 39.8, 42.9, 46.5, 45.9, 48.6, 48.2, 50.2, 50.8, 52.6  ] )     # Fitted Peak max values
TMax_Theory       = fTmax(Energy, E_c, X_0)                                                                                   # Get the TMax Variable from Theory
TMax_Ratio        = TMax_Theory / TMax_Truth_PFit # Residuals

# create a set of Parameters
p_tmax = Parameters()
p_tmax.add('X_0', value = X_0, min = 0,   max = 30,  vary = False)
p_tmax.add('E_c', value = E_c, min = 28,  max = 60,  vary = True)

# do fit, here with leastsq model
minner_tmax = Minimizer(fTmax_Residual, p_tmax, fcn_args=(Energy, TMax_Truth))
fitresult_tmax = minner_tmax.minimize()

# calculate the fit residuals
residual_tmax = fitresult_tmax.residual

# write error report
print "\n\n===== TMAX FIT REPORT ===== "
report_fit(fitresult_tmax)
print "===== ===================== ===== \n\n"
 
X_0_fit = fitresult_tmax.params['X_0'].value
E_c_fit = fitresult_tmax.params['E_c'].value 

# Generate the fit to plot 
TMax_Truth_fit = fTmax(Energy, E_c_fit,X_0_fit)

# Create a string for displaing on the plot
textstr_tmax = 'Fit Parameters \n---------------------\n $X_0=%.1f$ cm \n $ E_c=%.1f$ MeV \n $\chi_r^2=%.1f$ \n' % (X_0_fit, E_c_fit,fitresult_tmax.redchi )

# =+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+
# =+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+=+=++=+

# Moliere Radius
R_m_3D     = np.array([14.3, 13.7, 14.4, 13.9, 13.6, 13.7, 13.7, 13.7, 13.7, 13.7, 14.0, 13.7, 13.6, 13.5, 13.7])     # 3D Truth R_m
R_m_3D_2   = np.array([23.1, 22.2, 23, 22.7, 22.4, 22.4, 22.4, 22.8, 22.3, 22.5, 23.6, 22.6, 22.4, 22.3, 22.4 ])      # 3D Truth 2 * R_m
R_m_PDG    = 8.5  # PDG and Amaldi Definition of Moliere Radius
R_m_Fabjan = 11   # Fabjan and Amaldi Definition of Moliere Radius

# LArIAT in SBND Moliere Radius
R_m_LArIAT_SBND_3D = np.array([10.2, 9.5, 9.3, 8.8, 8.2, 8.2, 8, 7.9, 7.7, 7.7, 7.6, 7.4, 7.4, 7.2, 7.3])

R_m_LArIAT_3D = np.array([9.5, 9.7 , 9.3 , 9.09, 8.69]) # 3D truth moliere radius in LArIAT
Energy_LArIAT = np.array([100,200, 300, 500, 1000])  # Energy in LArIAT

# ----------------------- Plot the TMax Graph --------------------

plt.figure(0)
plt.plot(Energy , TMax_Truth_PFit , 'go',   label = 'Truth')
plt.plot(Energy, TMax_Truth_fit, 'b-', label = 'Fit')
plt.plot(Energy , TMax_Theory , 'r-',   label = 'Theory')
plt.text(2000, 15, textstr_tmax, fontsize = 15)
plt.ylabel("TMax [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() 
plt.savefig("./Plots/TMax.png", dpi=figqual)

plt.figure(1)
plt.plot(Energy , TMax_Ratio , 'bo',   label = 'TMax')
plt.ylabel("Truth/Theory [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() 
plt.savefig("./Plots/TMaxRatio.png", dpi=figqual)

# ----------------------- Plot the L @ 95% Graph --------------------
plt.figure(2)
plt.plot(Energy , Lon_Theory , 'r-',   label = 'Theory')
plt.plot(Energy , Lon_Truth , 'go',   label = 'Truth')
plt.plot(Energy, Lon_Truth_fit, 'b-', label = 'Fit')
plt.text(2200, 165, textstr_Lon, fontsize = 15)
plt.ylabel("Conatinment at 95% [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend() 
plt.savefig("./Plots/Lon_95%.png", dpi=figqual)

plt.figure(3)
plt.plot(Energy , Lon_Ratio , 'bo',   label = 'Lon(95%)')
plt.ylabel("Truth/Theory [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(fontsize = "small") 
plt.savefig("./Plots/Lon95%_Ratio.png", dpi=figqual)

# ----------------------- Section 2: Moliere Radius --------------------

plt.figure(4)
plt.plot(Energy , R_m_3D , 'ro',   label = 'Truth')
plt.axhline(y=R_m_Fabjan, color='b', linestyle='--', label = "Fabjan+ Amaldi")
plt.axhline(y=R_m_PDG, color='k', linestyle='--', label = "PDG + Amaldi")
plt.ylabel("Moliere Radius [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize = 13.5)
plt.savefig("./Plots/Moliere_Radius_SBND.png", dpi=figqual)

plt.figure(5)
plt.plot(Energy , 2 * R_m_3D , 'go',   label = 'Truth (2 * 90% Containment)')
plt.plot(Energy , R_m_3D_2 , 'ro',   label = 'Truth (95% Containment)')
plt.axhline(y=2 * R_m_Fabjan, color='b', linestyle='--', label = "Fabjan + Amaldi")
plt.axhline(y=2 * R_m_PDG, color='k', linestyle='--', label = "PDG + Amaldi")
plt.ylabel("2 * Moliere Radius [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize = 13.5)
plt.savefig("./Plots/Moliere_Radius_2_SBND.png", dpi=figqual)

# ----------------------- Plot the LArIAT Moliere Radius --------------------
plt.figure(6)
plt.plot(Energy ,R_m_LArIAT_SBND_3D , 'go',   label = 'Truth SBND(LArIAT geom)')
plt.plot(Energy_LArIAT , R_m_LArIAT_3D , 'ro',   label = 'Truth LArIAT')
plt.axhline(y= R_m_Fabjan, color='b', linestyle='--', label = "Fabjan + Amaldi")
plt.axhline(y=R_m_PDG, color='k', linestyle='--', label = "PDG + Amaldi")
plt.ylabel("Moliere Radius [cm]")
plt.xlabel("Energy [MeV]")
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0., fontsize = 13.5)
plt.savefig("./Plots/Moliere_Radius_LArIAT_SBND.png", dpi=figqual)