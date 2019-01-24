import numpy as np 
import matplotlib.pyplot as plt 

#data = np.loadtxt('gmin=122_tmax=3ns_11des2018/data_as_func_of_excitation_energy.dat.unchanged', skiprows=0, usecols=(0,1,2,3,4,5,6,7,8,9,10))
#unc = np.loadtxt('gmin=122_tmax=3ns_11des2018/uncertainties.dat', skiprows=0, usecols=(0,1,2))
data = np.loadtxt('gmin=122_tmax=3ns_noscale_14jan2019/data_as_func_of_excitation_energy.dat.unchanged', skiprows=0, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))
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
# plt.figure(0)
# #plt.errorbar(Ex, avg_ph_mult, yerr=unc[:,0], fmt="-x")
# plt.plot(Ex, avg_ph_mult_first, "ro-", label="First chance")
# plt.plot(Ex[3:-1], avg_ph_mult_second[3:-1], "ko-", label="Second chance")
# plt.plot(Ex, avg_ph_mult, "bo-", label="Total")
# #plt.title("Mg")
# #plt.title("Average photon multiplicity as function of Pu241 Ex-energy")
# plt.xlabel("Excitation energy [MeV]", fontsize=15)
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# plt.yticks(np.arange(6.9,7.4, step=0.05), [6.9,"",7.0,"",7.1,"",7.2,"",7.3,"",7.4], fontsize=14)
# plt.ylabel("Average photon multiplicity", fontsize=15)
# plt.legend()
# plt.grid()
# #plt.show()

# plt.figure(1)
# #plt.errorbar(Ex, avg_ph_energy, yerr=unc[:,1], fmt="-x")
# plt.plot(Ex, avg_ph_energy_first, "ro-", label="First chance")
# plt.plot(Ex[3:-1], avg_ph_energy_second[3:-1], "ko-", label="Second chance")
# plt.plot(Ex, avg_ph_energy, "bo-", label="Total")
# #plt.title("Average photon energy as function of Pu241 Ex-energy")
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# plt.yticks(np.arange(0.960,1, step=0.005), [0.960, "",0.970, "",0.980, "", 0.990, "", 1.00], fontsize=14)
# plt.xlabel("Excitation energy [MeV]", fontsize=15)
# plt.ylabel("Average photon energy [MeV]", fontsize=15)
# plt.legend()
# #plt.title("Eg")
# plt.grid()
# #plt.show()

# plt.figure(2)
# #plt.errorbar(Ex, total_ph_E, yerr=unc[:,2], fmt="-x")
# plt.plot(Ex, total_ph_E_first, "ro-", label="First chance")
# plt.plot(Ex[3:-1], total_ph_E_second[3:-1], "ko-", label="Second chance")
# plt.plot(Ex, total_ph_E, "bo-", label="Total")
# #plt.title("Etot")
# #plt.title("Total photon energy as function of Pu241 Ex-energy")
# plt.xlabel("Excitation energy [MeV]", fontsize=15)
# plt.ylabel("Total photon energy [MeV]", fontsize=15)
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# plt.yticks(np.arange(6.75, 7.2, step=0.05), ["", 6.8,"",6.9,"",7.0,"",7.1,"",7.2], fontsize=14)
# plt.legend()
# plt.grid()
# plt.show()

#With zoom-parts
f, axarr = plt.subplots(2, sharex=True)
plt.yticks(fontsize=14)
axarr[0].plot(Ex, avg_ph_energy_first, "ro-", label="First chance")
axarr[0].plot(Ex[3:-1], avg_ph_energy_second[3:-1], "ko-", label="Second chance")
axarr[0].plot(Ex, avg_ph_energy, "bo-", label="Total")
axarr[0].axis([5,13,0.6,1.3])
axarr[0].legend()
axarr[0].grid()
axarr[1].plot(Ex, avg_ph_energy_first, "ro-", label="First chance")
axarr[1].plot(Ex[3:-1], avg_ph_energy_second[3:-1], "ko-", label="Second chance")
axarr[1].plot(Ex, avg_ph_energy, "bo-", label="Total")
axarr[1].axis([5,13,0.955,1.002])
axarr[1].grid()
plt.yticks(fontsize=13)
plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
plt.xlabel("Excitation energy [MeV]", fontsize=15)
plt.ylabel("Average photon energy [MeV]", fontsize=15, position=(1,1))
plt.show()


# plt.plot(Ex, total_ph_E_first, "ro-", label="First chance", c= 'k')
# #plt.plot(Ex[3:-1], total_ph_E_second[3:-1], "ko-", label="Second chance")
# #plt.plot(Ex, total_ph_E, "bo-", label="Total")
# plt.axis([5.5,13,0,1.5])
# sub_axes = plt.axes([.6, .6, .25, .25])
# sub_axes.plot(Ex, total_ph_E_first, c = 'k')
# plt.grid()  
# plt.show()

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


# #Multichance fission handling
# first_second_data = np.genfromtxt('gmin=122_tmax=3ns_noscale_14jan2019/multichance_file.dat', skip_header=0, usecols=(0,1,2,3,4))

# Sn_Pu241 = 5.24152 #MeV
# Bf_240 = 6.05

# F = first_second_data[0][4]


# plt.plot(first_second_data[:,0] + Sn_Pu241, first_second_data[:,1]/F*100, "bx-", label="First chance")
# plt.plot(first_second_data[:,0] + Sn_Pu241, first_second_data[:,2]/F*100, "kx-", label="Second chance")
# plt.plot(first_second_data[:,0] + Sn_Pu241, first_second_data[:,3]/F*100, "gx-", label="Third chance")
# plt.axvline(x=Sn_Pu241 + Bf_240, color="r", label="$S_n(241Pu) + B_f(240Pu)$")
# plt.legend(fontsize=12)
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# plt.yticks(np.arange(0, 120, step=10), [0,"",20,"",40,"",60,"",80,"",100], fontsize=14)
# plt.ylabel("Share of fissions [%]", fontsize=14)
# plt.xlabel("Excitation energy [MeV]", fontsize=14)
# plt.grid()
# plt.show()

# #Fission fragment mass distribution
# fragment_mass_data = np.genfromtxt('gmin=122_tmax=3ns_noscale_14jan2019/fragment_mass_distr.dat', skip_header=0, usecols=(0,1,2))

# # plt.plot(fragment_mass_data[:,0] + Sn_Pu241, fragment_mass_data[:,1], label="FF1")
# # plt.plot(fragment_mass_data[:,0] + Sn_Pu241, fragment_mass_data[:,2], label="FF2")
# # plt.legend()
# # plt.xlabel("Ex= S_n + E_n [MeV]", fontsize=12)
# # plt.ylabel("Mass [A]")
# # plt.grid()
# # plt.show()

# f, (ax2, ax) = plt.subplots(2, 1, sharex=True)
# # plot the same data on both axes
# ax.plot(fragment_mass_data[:,0] + Sn_Pu241, fragment_mass_data[:,1], "bo-", label="Light fission fragment (as formed, not product)")
# ax2.plot(fragment_mass_data[:,0] + Sn_Pu241, fragment_mass_data[:,2], "ro-", label="Heavy fission fragment (as formed, not product)")
# ax.set_ylim(100, 102)  # outliers only
# ax2.set_ylim(139, 141)  # most of the data
# ax.spines['top'].set_visible(False)
# ax2.spines['bottom'].set_visible(False)
# plt.xlabel("Excitation energy [MeV]", fontsize=14)
# plt.ylabel("Mass [A]", position=(1,1))
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# ax.grid()
# ax.legend()
# ax2.grid()
# ax2.legend()
# plt.show()

# #Fission fragment kinetic energy
# fragment_kinE_data = np.genfromtxt('gmin=122_tmax=3ns_noscale_14jan2019/fragment_kinE.dat', skip_header=0, usecols=(0,1,2))
# f2, (ax22, ax2) = plt.subplots(2, 1, sharex=True)

# # plot the same data on both axes
# ax2.plot(fragment_kinE_data[:,0] + Sn_Pu241, fragment_kinE_data[:,1], "bo-", label="Light fission fragment")
# ax22.plot(fragment_kinE_data[:,0] + Sn_Pu241, fragment_kinE_data[:,2], "ro-", label="Heavy fission fragment")
# ax2.set_ylim(98, 102)  # outliers only
# ax22.set_ylim(72, 75)  # most of the data
# ax2.spines['top'].set_visible(False)
# ax22.spines['bottom'].set_visible(False)
# plt.xlabel("Excitation energy [MeV]", fontsize=14)
# plt.ylabel("Kinetic energy [MeV]", position=(1,1))
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# ax2.grid()
# ax2.legend()
# ax22.grid()
# ax22.legend()
# plt.show()

# #Neutron multiplicity

# neutron_data = np.genfromtxt('gmin=122_tmax=3ns_noscale_14jan2019/neutron_mult.dat', skip_header=0, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12))

# avg_n_mult = neutron_data[:,1]
# avg_n_energy = neutron_data[:,2]
# total_n_E = neutron_data[:,3]

# avg_n_mult_first = neutron_data[:,4]
# avg_n_energy_first = neutron_data[:,5]
# total_n_E_first = neutron_data[:,6]

# avg_n_mult_second = neutron_data[:,7]
# avg_n_energy_second = neutron_data[:,8]
# total_n_E_second = neutron_data[:,9]

# plt.plot(neutron_data[:,0] + Sn_Pu241, avg_n_mult_first, "ro-", label="First chance")
# plt.plot(neutron_data[:,0] + Sn_Pu241, avg_n_mult_second, "ko-", label="Second chance")
# plt.plot(neutron_data[:,0] + Sn_Pu241, avg_n_mult, "bo-", label="Total")
# plt.xlabel("Excitation energy [MeV]", fontsize=15)
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# #plt.yticks(np.arange(6.9,7.4, step=0.05), [6.9,"",7.0,"",7.1,"",7.2,"",7.3,"",7.4], fontsize=14)
# plt.ylabel("Average neutron multiplicity", fontsize=15)
# plt.legend()
# plt.grid()
# plt.show()

# plt.plot(neutron_data[:,0] + Sn_Pu241, total_n_E_first, "ro-", label="First chance")
# plt.plot(neutron_data[:,0] + Sn_Pu241, total_n_E_second, "ko-", label="Second chance")
# plt.plot(neutron_data[:,0] + Sn_Pu241, total_n_E, "bo-", label="Total")
# plt.xlabel("Excitation energy [MeV]", fontsize=15)
# plt.ylabel("Total neutron energy [MeV]", fontsize=15)
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# #plt.yticks(np.arange(6.75, 7.2, step=0.05), ["", 6.8,"",6.9,"",7.0,"",7.1,"",7.2], fontsize=14)
# plt.legend()
# plt.grid()
# plt.show()

# plt.plot(neutron_data[:,0] + Sn_Pu241, avg_n_energy_first, "ro-", label="First chance")
# plt.plot(neutron_data[:,0] + Sn_Pu241, avg_n_energy_second, "ko-", label="Second chance")
# plt.plot(neutron_data[:,0] + Sn_Pu241, avg_n_energy, "bo-", label="Total")
# plt.xticks(np.arange(5, 13, step=0.5), [5," " ,6, " ",7, " ",8, " ",9, " " ,10, " ",11, " ",12, " ",13], fontsize=14)
# #plt.yticks(np.arange(0.960,1, step=0.005), [0.960, "",0.970, "",0.980, "", 0.990, "", 1.00], fontsize=14)
# plt.xlabel("Excitation energy [MeV]", fontsize=15)
# plt.ylabel("Average neutron energy [MeV]", fontsize=15)
# plt.legend()
# #plt.title("Eg")
# plt.grid()
# plt.show()


