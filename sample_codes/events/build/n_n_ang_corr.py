#Calculates angular correlation between emitted n's
#must first have run EventToRoot_compilable.C, as the input file ndir.txt is made there

import numpy as np 
import matplotlib.pyplot as plt

infile = open('n_dir.txt', 'r')

cos_angles = []

for line in infile:
	words = line.split()
	event_number = float(words[0])
	list1 = []
	for k in range(1, len(words)):
		list1.append(float(words[k]))
	
	nu = len(list1)/3
	p_xyz = np.zeros((nu,3))
	#want to put all neutrons in nested array p_xyz[n][0:2] is the x,y,z-direction of neutron n

	for i in range(0,nu):
		p_xyz[i][0] = list1[3*i]
		p_xyz[i][1] = list1[3*i+1]
		p_xyz[i][2] = list1[3*i+2]

	for i in range(nu-1):
		for j in range(i+1,nu,1):
			total = np.dot(p_xyz[i][0:3],p_xyz[j][0:3])/(np.sqrt(p_xyz[i][0]**2+p_xyz[i][1]**2+p_xyz[i][2]**2)*np.sqrt(p_xyz[j][0]**2+p_xyz[j][1]**2+p_xyz[j][2]**2))
			cos_angles.append(total)
			


plt.hist(cos_angles, bins=100)
plt.ylabel('Number of events')
plt.xlabel('cosine of angle')
plt.show()
