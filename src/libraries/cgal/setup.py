# python setup.py build build_ext --inplace
# python setup.py build
# python steup.py install

from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
    #ext_modules=[Extension("_chi2", ["_chi2.cpp", "chi.cpp"])],
    ext_modules=[
        Extension("elipsoid", ["_elipsoid.cpp"],
                  library_dirs = ["."],
                  libraries = ["elipsoid"],
                  #runtime_library_dirs=["$ORIGIN/cgal"]
                       )],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)

