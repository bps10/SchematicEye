#!/usr/bin/env python
import sys
sys.path.append("./")
import cProfile
import timeit
import eye as eye



print (timeit.timeit('eye.py_eye( object_distance=1e8, off_axis=0, pupil_size=4, diopters=0)', 
	number=10, setup='import eye as eye'))

intensity = eye.py_eye( object_distance=1e8, off_axis=0, 
		pupil_size=4, diopters=0)

import matplotlib.pylab as plt
plt.figure()
plt.plot(intensity)
plt.show()