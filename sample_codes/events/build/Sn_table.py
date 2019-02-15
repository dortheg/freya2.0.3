import matplotlib.pyplot as plt 
import numpy as np

Sn_table = np.zeros((119,177))

Z = np.genfromtxt("../../../data_freya/MassAudi.dat", skip_header=5, usecols=(0))
A = np.genfromtxt("../../../data_freya/MassAudi.dat", skip_header=5, usecols=(1))
B = np.genfromtxt("../../../data_freya/MassAudi.dat", skip_header=5, usecols=(4))

#print(Sn_table[118][176])

for i in range(1,len(Z),1):
	if Z[i] == Z[i-1]:
		S_n = B[i] - B[i-1]
		z = int(Z[i])
		n = int(A[i] - Z[i])
		Sn_table[z][n] = S_n

