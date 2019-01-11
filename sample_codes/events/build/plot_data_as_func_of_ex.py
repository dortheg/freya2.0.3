import numpy as np 
import matplotlib.pyplot as plt 

#data = np.loadtxt('gmin=122_tmax=3ns_11des2018/data_as_func_of_excitation_energy.dat.unchanged', skiprows=0, usecols=(0,1,2,3,4,5,6,7,8,9,10))
#unc = np.loadtxt('gmin=122_tmax=3ns_11des2018/uncertainties.dat', skiprows=0, usecols=(0,1,2))
data = np.loadtxt('data_as_func_of_excitation_energy.dat', skiprows=0, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))
#n0 = np.loadtxt('pre-fission-neutrons.dat', skiprows=0, usecols=(0,1,2))

Ex = data[:,4]
#Ex_1 = [0,1,2,3,4,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,8.0]

avg_ph_mult = data[:,5]
avg_ph_energy = data[:,7]
total_ph_E = data[:,9]

avg_ph_mult_first = data[:,11]
avg_ph_energy_first = data[:,13]
total_ph_E_first = data[:,15]

avg_ph_mult_second = data[:,17]
avg_ph_energy_second = data[:,19]
total_ph_E_second = data[:,21]

#avg_ph_mult_third = data[:,23]
#avg_ph_energy_third = data[:,25]
#total_ph_E_third = data[:,27]

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#plt.plot(Ex, avg_ph_mult, 'o')
plt.figure(0)
#plt.errorbar(Ex, avg_ph_mult, yerr=unc[:,0], fmt="-x")
plt.plot(Ex, avg_ph_mult, "bo-", label="Total")
plt.plot(Ex, avg_ph_mult_first, "ro-", label="First")
plt.plot(Ex, avg_ph_mult_second, "ko-", label="Second")
#plt.title("Average photon multiplicity as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy [MeV]")
plt.ylabel("Average photon multiplicity")
plt.legend()
plt.grid()
#plt.show()

plt.figure(1)
#plt.errorbar(Ex, avg_ph_energy, yerr=unc[:,1], fmt="-x")
plt.plot(Ex, avg_ph_energy, "bo-", label="Total")
plt.plot(Ex, avg_ph_energy_first, "ro-", label="First")
plt.plot(Ex, avg_ph_energy_second, "ko-", label="Second")
#plt.title("Average photon energy as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy [MeV]")
plt.ylabel("Average photon energy [MeV]")
plt.legend()
plt.grid()
#plt.show()

plt.figure(2)
#plt.errorbar(Ex, total_ph_E, yerr=unc[:,2], fmt="-x")
plt.plot(Ex, total_ph_E, "bo-", label="Total")
plt.plot(Ex, total_ph_E_first, "ro-", label="First")
plt.plot(Ex, total_ph_E_second, "ko-", label="Second")
#plt.title("Total photon energy as function of Pu241 Ex-energy")
plt.xlabel("Excitation energy [MeV]")
plt.ylabel("Total photon energy [MeV]")
plt.legend()
plt.grid()
plt.show()

# plt.figure(3)
# plt.plot(Ex_1, n0[:,0], '-o')
# plt.plot(Ex_1, n0[:,1], '-o')
# plt.plot(Ex_1, n0[:,2], '-o')
# plt.title("Avg number of pre-fission neutrons as function of Pu241 Ex-energy")
# plt.axis([0,max(Ex),0,2.5])
# #plt.xticks(np.arange(min(Ex), max(Ex), 0.5)) #makes ticks alog x-axis
# #plt.yticks(np.arange(min(n0[:,0]), max(n0[:,0]), 0.1)) #ticks along y-axis
# plt.grid(True)
# plt.xlabel("Excitation energy of Pu241*")
# plt.ylabel("Avg number of pre-fission neutrons")
# plt.show()

