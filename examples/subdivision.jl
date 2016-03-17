using SurfaceGeometry

isdefined(:ViewerActive) || (ViewerActive=false)

# using JLD
# @load "sphere.jld"

fdis(x,y,z) = x^2 + y^2 + z^2 - 1
points, faces = SurfaceMesh(fdis,CGALSurfaceMesher())


#rpoints,rfaces = subdivision(points,faces;method=:linear)
#rpoints,rfaces = subdivision(points,faces;method=:paraboloid)
rpoints,rfaces = subdivision(points,faces,x -> x[1]^2 + x[2]^2 + x[3]^2 - 1)

points,faces  = rpoints, rfaces

using Escher
if !ViewerActive
    include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
    @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
    ViewerActive = true
end




