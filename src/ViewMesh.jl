#!/usr/bin/env julia
# It would make sense to run it in REPL so no need to fnames would rise

# module ViewMesh
# export view3d

function ConvertToTuple(p,t)
    points = Array(Tuple{Float64,Float64,Float64},size(p,2))
    for i in 1:size(p,2)
        points[i] = tuple(p[:,i]...)
    end
    
    faces = Array(Tuple{Int64,Int64,Int64},size(t,2))
    for i in 1:size(t,2)
        faces[i] = tuple(t[:,i]...)
    end

    return points,faces
end

# The function which should be run on another process
#using JLD

#points, faces = loadmesh

using Escher

function view3d()

    # eval(:(using Escher))
    # eval(:(import ThreeJS))
    # eval(:(using ThreeJS))
    # eval(:(using Compat))
    
    # For now a simple test
    
    include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
    escher_serve(5555,Pkg.dir("SurfaceGeometry","src","libraries","viewers"))
end

# end

# using ViewMesh
