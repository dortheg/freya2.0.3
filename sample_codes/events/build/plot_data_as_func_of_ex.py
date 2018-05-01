import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('data_as_func_of_excitation_energy.dat.unchanged', skiprows=0, usecols=(0,1,2,3,4,5,6,7,8,9,10))
#unc = np.loadtxt('uncertainties.dat', skiprows=0, usecols=(0,1,2))
n0 = np.loadtxt('pre-fission-neutrons.dat', skiprows=0, usecols=(0,1,2))

#Ex = data[:,4]
Ex = [0,1,2,3,4,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,8.0]
"""
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
"""
plt.figure(3)
plt.plot(Ex, n0[:,0], '-o')
plt.plot(Ex, n0[:,1], '-o')
plt.plot(Ex, n0[:,2], '-o')
plt.title("Avg number of pre-fission neutrons as function of Pu241 Ex-energy")
plt.axis([0,max(Ex),0,2.5])
#plt.xticks(np.arange(min(Ex), max(Ex), 0.5)) #makes ticks alog x-axis
#plt.yticks(np.arange(min(n0[:,0]), max(n0[:,0]), 0.1)) #ticks along y-axis
plt.grid(True)
plt.xlabel("Excitation energy of Pu241*")
plt.ylabel("Avg number of pre-fission neutrons")
plt.show()

