# CFLAGS=-I$(IJULIA) -I$(ICxxWrap) -L$(LDIR)
#eltopo: eltopo.cpp
#	g++ eltopo.cpp -fPIC -std=c++11 -shared -o eltopo.so

eltopo: 
	g++ eltopo.cpp libeltopo_release.a -o eltopo.so -fPIC -std=c++11 -shared -llapack -lblas -lstdc++ -lm -I../common -I../eltopo3d 
