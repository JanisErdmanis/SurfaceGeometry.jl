using SurfaceGeometry
using Base.Test

# Pkg.test("SurfaceGeometry")
# write your own tests here
#@test 1 == 2

### Loading of spherical mesh

function SphereError(points)
    error = 0
    for xkey in 1:size(points,2)
        pos = points[:,xkey]
        n = pos/norm(pos)
        err = norm(n - pos)
        if err>error
            error = err
        end
    end
    return error
end

### Loading of mesh for other tests
println("Mesh loading test")
using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

@test SphereError(points) < 0.01

println("Mesh generation with CGAL")

### Interface will change to pass signed distance function
### CGAL mesh generator
points, faces = ellipsoid_mesh_cgal(1,1,1)
@test SphereError(points) < 0.01
    
### Matlab mesh generator
println("Mesh generation with distmesh")
try
    points, faces = ellipsoid_mesh_matlab(1,1,1,0.35)
    @test SphereError(points) < 0.01
catch
    warn("does not work")
end

### Testing a pushback
sdist(x) = x[1]^2 + x[2]^2 + x[3]^2 - 1
@test isapprox(norm(pushback(sdist,[1.,1.,1.])),1)

### Testing topology

println("Topology function tests")

t = [5,7,10,7,5,6,4,0,3,0,4,6,4,7,6,4,9,10,7,4,10,0,2,1,2,0,6,2,5,1,5,2,6,8,4,3,4,11,9,8,11,4,9,11,3,11,8,3]
faces = reshape(t,(3,div(length(t),3))) + 1 

triangles = []
for i in FaceRing(1,faces)
    #println("i")
    push!(triangles,i)
end

@test sort(triangles)==[3,4,8,9]

##### Surface Properties ######

println("Tests for surface properties")

using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

angle = 0
curvaturer = 0

for xkey in 1:size(points,2)
    iter = FaceRing(xkey,faces)
    ncalc = vnormal(xkey,points,faces,iter)
    ccalc = vcurvature(xkey,points,faces,iter,ncalc)

    nT = points[:,xkey]/norm(points[:,xkey])
    angl = acos(dot(ncalc,nT))
    if angl>angle
        angle = angl
    end

    cT = 1.
    N = size(points,2)
    curvaturer += abs((ccalc - cT)/cT)/N
end

@test angle*180/pi < 1.3
@test curvaturer < 0.01

@test isapprox(volume(points,faces),4/3*pi,atol=0.05)
@test isapprox(*(ellipsoid_parameters(points)...),1)

### And now integrators

println("With Eigene test checking surface integrators")

N = 100
meth = AdamsStep
PasiveStabilisation = false
include("integrator.jl")
@test SphereError(points) < 0.0001

N = 10000
meth = EilerStep
PasiveStabilisation = false
include("integrator.jl")
@test SphereError(points) < 0.001

N = 100
meth = AdamsStep
PasiveStabilisation = true
include("integrator.jl")
@test SphereError(points) < 0.0005

