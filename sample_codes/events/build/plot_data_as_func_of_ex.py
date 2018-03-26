import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('data_as_func_of_excitation_energy.dat.unchanged', skiprows=0, usecols=(0,1,2,3,4,5,6,7,8,9,10))
unc = np.loadtxt('uncertainties.dat', skiprows=0, usecols=(0,1,2))

Ex = data[:,4]
avg_ph_mult = data[:,5]
avg_ph_energy = data[:,7]
total_ph_E = data[:,9]

#plt.plot(Ex, avg_ph_mult, 'o')
plt.figure(0)
plt.errorbar(Ex, avg_ph_mult, yerr=unc[:,0], fmt="x")
#plt.plot(Ex, avg_ph_mult, "o")
plt.title("Average photon multiplicity as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy of Pu241*")
plt.ylabel("Average photon multiplicity")
#plt.show()

plt.figure(1)
plt.errorbar(Ex, avg_ph_energy, yerr=unc[:,1], fmt="x")
#plt.plot(Ex, avg_ph_energy, "o")
plt.title("Average photon energy as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy of Pu241*")
plt.ylabel("Average photon energy")
#plt.show()

plt.figure(2)
#plt.errorbar(Ex, total_ph_E, yerr=unc[:,2], fmt="x")
plt.plot(Ex, total_ph_E, "o")
plt.title("Total photon energy as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy of Pu241*")
plt.ylabel("Average photon energy")
plt.show()

