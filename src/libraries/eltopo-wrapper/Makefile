ElTopoDIR = ../eltopo-master

copydefs:
	cp Makefile.local_defs $(ElTopoDIR)/eltopo3d/

eltopo: copydefs
	cd $(ElTopoDIR)/eltopo3d && make depend && make release

eltopowrap: eltopo.cpp eltopo
	g++ eltopo.cpp $(ElTopoDIR)/eltopo3d/libeltopo_release.a -o eltopo.so -fPIC -shared -llapack -lblas -lstdc++ -lm -I$(ElTopoDIR)/common -I$(ElTopoDIR)/eltopo3d 

