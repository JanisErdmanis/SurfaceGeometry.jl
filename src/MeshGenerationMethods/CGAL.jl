# include("Utils.jl")

immutable CGALSurfaceMesher
    AngularBound::Float64
    RadiusBound::Float64
    DistanceBound::Float64
    BoundingRadius::Float64
end

### For an easy syntax
function CGALSurfaceMesher(;AngularBound=20,RadiusBound=0.5,DistanceBound=0.1,BoundingRadius=8.0)
    CGALSurfaceMesher(AngularBound,RadiusBound,DistanceBound,BoundingRadius)
end

# I will try to play around with Cxx otherwise I can directly interface Cxx code
# Some kind of mesh simplificaction is needed for using it in computations
function EllipsoidMesh(a,b,c,cg::CGALSurfaceMesher)

    ENV["PYTHON"] = "/usr/bin/python3" 
    eval(:(Pkg.build("PyCall")))
    eval(:(using PyCall))
    unshift!(PyVector(pyimport("sys")["path"]), "")

    olddir = pwd()
    try 
        cd(Pkg.dir("SurfaceGeometry","src","libraries"))
        eval(:(@pyimport cgal))
    finally
        cd(olddir)
    end
    
    p, t = cgal.elipsoid_mesh(a,b,c,cg.AngularBound,cg.RadiusBound,cg.DistanceBound,cg.BoundingRadius)
    p = transpose(p)
    t = transpose(t+1) # since fortran indexing

    if !isoriented(p,t)
        t[1,:], t[2,:] = t[2,:], t[1,:]
    end

    p,t = map(Float64,p), map(Int,t)

    for vi in 1:size(p)[end]
        p[:,vi] = pushback(x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1,p[:,vi])
    end
    
    return p,t
end

SurfaceMesh(fdis::Function,ds::CGALSurfaceMesher) = error("Not implemented")
