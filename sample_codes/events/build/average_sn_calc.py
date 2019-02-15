import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

Z1 = np.genfromtxt("Product_nuclei.dat", usecols=(0))
A1 = np.genfromtxt("Product_nuclei.dat", usecols=(1))

Z2 = np.genfromtxt("Product_nuclei.dat", usecols=(2))
A2 = np.genfromtxt("Product_nuclei.dat", usecols=(3))

#Create binding energy table
Sn_table = np.zeros((119,177))

Z_tab = np.genfromtxt("../../../data_freya/MassAudi.dat", skip_header=5, usecols=(0))
A_tab = np.genfromtxt("../../../data_freya/MassAudi.dat", skip_header=5, usecols=(1))
B_tab = np.genfromtxt("../../../data_freya/MassAudi.dat", skip_header=5, usecols=(4))


for i in range(1,len(Z_tab),1):
	if Z_tab[i] == Z_tab[i-1]:
		S_n = B_tab[i] - B_tab[i-1]
		z = int(Z_tab[i])
		n = int(A_tab[i] - Z_tab[i])
		Sn_table[z][n] = S_n

fissprod = np.zeros((94,146))

for i in range(len(Z1)):
	z1 = int(Z1[i])
	n1 = int(A1[i] - Z1[i])

	z2 = int(Z2[i])
	n2 = int(A2[i] - Z2[i])

	fissprod[z1][n1] += 1
	fissprod[z2][n2] += 1

#Calculate average separation energy
S_n_tot = 0
counter = 0
for i in range(94):
	for j in range(146):
 		if fissprod[i][j]>0:
 			S_n_tot += fissprod[i][j]*Sn_table[i][j]
 			#print(i,j,Sn_table[i][j])
 			counter += fissprod[i][j]


S_n_tot = S_n_tot/(2*len(Z1)) #Because two fission fragments
print(S_n_tot)


#Plot fission fragment distr
#plt.imshow(fissprod)
#plt.show()


