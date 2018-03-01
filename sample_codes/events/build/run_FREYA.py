import os

"""
Program that runs freya, creates Pu240.dat FREYA-file, then runs EventToRoot_compilable and creates Pu240.dat.root, 
which is finally analyzed with root.

All interesting quantities are printed in the file data_as_func_of_excitation_enery.dat

OBS:: must run: g++ EventToRoot_compilable.C -o EventToRoot_compilable `root-config --cflags --libs` -std=c++14 when 
changing the filename in EventToRoot_compilable.C, as this program is unable to do that by itself

OBS:: Must remember to change number of fissions in freya_root_analyzer.C, if else that 100k! Else uncertainties are wrong
"""

from subprocess import Popen, PIPE

import numpy as np 

Ex = np.linspace(0,6,10)

for i in range(len(Ex)):
	p = Popen('./events', stdin=PIPE)
	p.communicate(os.linesep.join(["94", "240", "1", "%f" % Ex[i], "10000", "Pu240.dat"]))

	o = Popen("./EventToRoot_compilable", stdin=PIPE)
	o.communicate(os.linesep.join([" "])) #makes sure the program stops?

	#m = Popen('root', stdin=PIPE)
	#m.communicate(os.linesep.join(['.x freya_root_analyzer.C']))

	from subprocess import call
	call(["root","-q", "-l","freya_root_analyzer.C"])
	#call(["root","-q", "-l","fission_script_root6.C"])


