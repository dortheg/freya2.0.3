#This script reads avg and total gamma energy and gamma multiplicity from the file avg_values_for_unc_cal.dat
#Calculates the statistical uncertainties in these values for one excitation energy -> must be done for all excitation energies

import numpy as np
import matplotlib.pyplot as plt 

data = np.loadtxt('avg_values_for_unc_calc.dat', skiprows=0, usecols=(0,1,2))

#Photon multiplicity
N_g = data[:,0]

#Avg gamma energy
E_g = data[:,1]

#Total gamma energy
E_g_tot = data[:,2]

#Create arrays
N_g = np.array(N_g)
E_g = np.array(E_g)
E_g_tot = np.array(E_g_tot)

#Uncertainty calc
sigma_N_g = 0
sigma_E_g = 0
sigma_E_g_tot = 0

for i in range(0, len(N_g)):
	sigma_N_g += (N_g[i] - np.average(N_g))**2
	sigma_E_g += (E_g[i] - np.average(E_g))**2
	sigma_E_g_tot += (E_g_tot[i] - np.average(E_g_tot))**2

sigma_N_g = np.sqrt(sigma_N_g/(len(N_g) - 1))
sigma_E_g = np.sqrt(sigma_E_g/(len(E_g) - 1))
sigma_E_g_tot = np.sqrt(sigma_E_g_tot/(len(E_g_tot) - 1))

infile = open('uncertainties.dat', "a")
infile.write("%f  %f  %f \n" % (sigma_N_g, sigma_E_g, sigma_E_g_tot))


