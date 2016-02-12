cd(Pkg.dir("SurfaceGeometry","src","libraries","cgal"))
run(`g++ -fPIC -shared elipsoid.cpp -o libelipsoid.so -lmpfr -lgmp -lCGAL -lboost_thread`)

