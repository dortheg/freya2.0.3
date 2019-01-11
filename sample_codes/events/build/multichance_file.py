import matplotlib.pyplot as plt 
import numpy as np 

first_second_data = np.genfromtxt('multichance_file.dat', skip_header=0, usecols=(0,1,2,3,4))

Sn_Pu241 = 5.24152 #MeV
Bf_240 = 6.05

F = first_second_data[0][4]


plt.plot(first_second_data[:,0] + Sn_Pu241, first_second_data[:,1]/100, label="First chance")
plt.plot(first_second_data[:,0] + Sn_Pu241, first_second_data[:,2]/100, label="Second chance")
plt.plot(first_second_data[:,0] + Sn_Pu241, first_second_data[:,3]/100, label="Third chance")
plt.axvline(x=Sn_Pu241 + Bf_240, color="r", label="Sn(241Pu) + Bf(240Pu)")
plt.legend(fontsize=12)
plt.ylabel("Percent of fissions", fontsize=12)
plt.xlabel("Ex= S_n + E_n [MeV]", fontsize=12)
plt.grid()
plt.show()
