#!/usr/bin/env python
import distutils
from distutils.core import setup

description = "Compute likelihood given data, model and emulator"

setup(name="lya_likelihood", 
    version="0.1.0",
    description=description,
    url="https://github.com/andreufont/LyaCosmoParams/tree/master/lya_likelihood",
    author="Andreu Font-Ribera",
    py_modules=['lya_theory','compute_likelihood'],
    package_dir={'': 'py'})

