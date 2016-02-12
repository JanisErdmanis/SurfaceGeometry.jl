all: pymodule

# export LD_LIBRARY_PATH="./cgal"

CXX_FLAGS = -frounding-math
CXX_DEFINES = -DCGAL_USE_GMP -DCGAL_USE_MPFR -DCGAL_USE_ZLIB
CXX_LIBS = -lmpfr -lgmp -lCGAL -lboost_thread

# # %.o: $(wildcard *.cpp)
# # 	c++ $(CXX_DEFINES) $(CXX_FLAGS) -o $@ -c $^

### C++ libary

libelipsoid.so: elipsoid.cpp
	g++ -fPIC -shared $^ -o $@ -lmpfr -lgmp -lCGAL -lboost_thread #-std=c++11

# surface: testing.cpp libelipsoid.so
# 	c++ testing.cpp libelipsoid.so -o $@ #-L. -l:elipsoid.so

### Wrapping by hand

pymodule: _elipsoid.cpp libelipsoid.so setup.py
	LDFLAGS="-Wl,-rpath ." \
	python setup.py build build_ext --inplace


### SWIG

TCL_INCLUDE   = -I/usr/local/include
TCL_LIB       = -L/usr/local/lib

elipsoid_wrap.cxx: elipsoid.cpp elipsoid.i
	swig -c++ -tcl8 elipsoid.i	

elipsoid.o: elipsoid.cpp
	c++ -fpic $(CXX_DEFINES) $(CXX_FLAGS) -o $@ -c $^ -I/usr/include/python2.7

elipsoid_wrap.o: elipsoid_wrap.cxx
	c++ -fpic $(CXX_DEFINES) $(CXX_FLAGS) -o $@ -c $^ -I/usr/include/tcl8.4 -I/usr/include/python2.7

elipsoid_swig.so: elipsoid.o elipsoid_wrap.o
	c++ -fpic -shared $^ -o $@ -lmpfr -lgmp -lCGAL -lboost_thread


## python bindings from swig

elipsoid.py: elipsoid.cpp elipsoid_wrap.cxx elipsoid.i
	swig -python elipsoid.i

_elipsoid.so: elipsoid.o elipsoid_wrap.o
	ld -shared $^ -o $@ #-lmpfr -lgmp -lCGAL -lboost_thread


## cleaning
clean:
	rm -rf elipsoid.cpython-33m.so;
	rm -rf build


