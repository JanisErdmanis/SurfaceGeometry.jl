cd(joinpath(dirname(dirname(@__FILE__)),"src","libraries","cgal"))
run(`g++ -fPIC -shared cgal_mesh.cpp -o libcgalmesh.so -lmpfr -lgmp -lCGAL -lboost_thread -frounding-math`)
cd(joinpath(dirname(dirname(@__FILE__)),"src","libraries","eltopo-wrapper"))
run(`make eltopowrap`)
