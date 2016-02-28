cd(Pkg.dir("SurfaceGeometry","src","libraries","cgal"))
run(`g++ -fPIC -shared cgal_mesh.cpp -o libcgalmesh.so -lmpfr -lgmp -lCGAL -lboost_thread -frounding-math`)
