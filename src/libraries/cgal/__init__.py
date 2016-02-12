import os

initdir = os.path.dirname(os.path.realpath(__file__))
curdir = os.getcwd()

os.chdir(initdir)
os.system("LDFLAGS=\"-Wl,-rpath .\" python setup.py build build_ext --inplace")
from elipsoid import *

os.chdir(curdir)
