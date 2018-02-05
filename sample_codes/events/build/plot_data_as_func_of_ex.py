import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('data_as_func_of_excitation_energy.dat', skiprows=6, usecols=(0,1,2,3,4,5,6,7,8,9,10))

Ex = data[:,4]
avg_ph_mult = data[:,5]
avg_ph_energy = data[:,7]

#plt.plot(Ex, avg_ph_mult, 'o')
plt.errorbar(Ex, avg_ph_mult, yerr=data[:,6], fmt="o")
plt.title("Average photon multiplicity as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy of Pu241*")
plt.ylabel("Average photon multiplicity")
plt.show()

plt.errorbar(Ex, avg_ph_energy, yerr=data[:,8], fmt="o")
plt.title("Average photon energy as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy of Pu241*")
plt.ylabel("Average photon energy")
plt.show()