### This file is meant for small functions and to interface more complicated code


using ForwardDiff               # Can be made optional

function pushback(sdist::Function,x::Array{Float64,1})
    for i in 1:100
        grad = ForwardDiff.gradient(sdist,x)
        x = x - sdist(x)*grad/norm(grad)^2
        # I allow the root to be only on gradient
    end
    return x
end    


### Loading mesh by filename in the /preproc directgory
### At this point it should be unimportant if it is either generated with cgal or matlab


### Creation routines


function isoriented{T <: AbstractFloat, S <: Integer}(points::Array{T,2},triangles::Array{S,2})

    tri = triangles[:,1]
    p1 = points[:,tri[1]]
    p2 = points[:,tri[2]]
    p3 = points[:,tri[3]]

    a1 = p2 - p1
    a2 = p3 - p1

    normal = cross(a1,a2)
    normal /= sqrt(dot(normal,normal))

    direction = p1
    direction /= sqrt(dot(p1,p1))

    angle = acos(dot(normal,direction))

    if angle<pi/2
        return true
    else
        return false
    end
end

# should be rewritten in Julia at some point 
function ellipsoid_mesh_matlab(a::Real,b::Real,c::Real,step::Real)

    olddir = pwd()
    try
        cd(Pkg.dir("SurfaceGeometry","src","libraries","distmesh"))
        Base.syntax_deprecation_warnings(false)
        eval(:(using MATLAB))

        a = float(a)
        b = float(b)
        c = float(c)
        step = float(step)
        p, t = mxcall(:elipsoid,2,a,b,c,step)
    finally
        cd(olddir)
    end
        
    
    t = transpose(map(Int,t))
    p = transpose(p)

    ### Converting to tuples
    ### Checking if there is need to change order for triangles    


    if !isoriented(p,t)
        t[1,:], t[2,:] = t[2,:], t[1,:]
    end

    p,t = map(Float64,p), map(Int,t)
    
    for vi in 1:size(p)[end]
        p[:,vi] = pushback(x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1,p[:,vi])
    end
    
    return p,t
end

# I will try to play around with Cxx otherwise I can directly interface Cxx code
# Some kind of mesh simplificaction is needed for using it in computations
function ellipsoid_mesh_cgal(a,b,c,angular_bound=20,radius_bound=0.5,distance_bound=0.1,bounding_radius=8.0)

    eval(:(using PyCall))
    unshift!(PyVector(pyimport("sys")["path"]), "")

    olddir = pwd()
    try 
        cd(Pkg.dir("SurfaceGeometry","src","libraries"))
        eval(:(@pyimport cgal))
    finally
        cd(olddir)
    end
    
    p, t = cgal.elipsoid_mesh(a,b,c,angular_bound,radius_bound,distance_bound,bounding_radius)
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

# Making surface mesh directly from octahedron
# Subdivision scheme which should be easy to implement as I already have pushback


