#!/usr/local/bin/python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np 

setup(
	ext_modules=[Extension("eye", ["python/eye.pyx", "src/eye_eye.cpp", "src/eye_analysis.cpp"], 
	libraries=["Goptical"],
	language="c++",)],
	include_dirs = [np.get_include()],
	cmdclass = {'build_ext': build_ext})
