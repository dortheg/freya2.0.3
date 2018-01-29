import os

"""
Program that runs freya, creates Pu240.dat FREYA-file, then runs EventToRoot_compilable and creates Pu240.dat.root, 
which is finally analyzed with root.

All interesting 
"""

from subprocess import Popen, PIPE



for i in range(1):
	p = Popen('./events', stdin=PIPE)
	p.communicate(os.linesep.join(["94", "240", "1", "5", "1000", "Pu240.dat"]))

	o = Popen("./EventToRoot_compilable", stdin=PIPE)
	o.communicate(os.linesep.join([" "])) #makes sure the program stops?

	#m = Popen('root', stdin=PIPE)
	#m.communicate(os.linesep.join(['.x freya_root_analyzer.C']))

	from subprocess import call
	call(["root","-q","freya_root_analyzer.C"])

