using SurfaceGeometry
using LinearAlgebra
using Test

# Pkg.test("SurfaceGeometry")
# write your own tests here
# @test 1 == 2

# Loading of spherical mesh

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
@info "Mesh loading test"
using JLD
data = load(joinpath(dirname(@__FILE__),"sphere.jld"))
points = data["points"]
faces = data["faces"]

@test SphereError(points) < 0.01

@info "Topology function tests"

t = [5,7,10,7,5,6,4,0,3,0,4,6,4,7,6,4,9,10,7,4,10,0,2,1,2,0,6,2,5,1,5,2,6,8,4,3,4,11,9,8,11,4,9,11,3,11,8,3]
faces = reshape(t,(3,div(length(t),3))) .+ 1 

triangles = []
for i in FaceVRing(5,faces)
    #println("i")
    push!(triangles,i)
end

@test sort(triangles)==[3,4,5,6,7,12,13,14]

triverticies = []
for i in DoubleVertexVRing(5,faces)
    #println("i")
    push!(triverticies,i)
end

@test (1,4) in triverticies

verticies = []
for i in VertexVRing(5,faces)
    push!(verticies,i)
end

@test sort(verticies)==[1,4,7,8,9,10,11,12] #[1,2,6,7]

@info "Testing Complex DS"

points = zero(faces)
fb = FaceBasedDS(faces)

triangles = []
for i in FaceVRing(5,fb)
    push!(triangles,i)
end
@test sort(triangles)==[3,4,5,6,7,12,13,14]

triverticies = []
for i in DoubleVertexVRing(5,fb)
    push!(triverticies,i)
end

@test (1,4) in triverticies

verticies = []
for i in VertexVRing(5,fb)
     push!(verticies,i)
end

@test sort(verticies)==[1,4,7,8,9,10,11,12]

@info "Topology tests for Connectivity table"

using JLD
data = load(joinpath(dirname(@__FILE__),"sphere.jld"))
#data = load("sphere.jld")
faces = data["faces"]

con = ConnectivityDS(faces,10)

verticies = []
for i in VertexVRing(1,con)
    push!(verticies,i)
end

@test sort(verticies)==[2,4,5,14,15]

triverticies = []
for i in DoubleVertexVRing(1,con)
    push!(triverticies,i)
end

@test (2,14) in triverticies

# Here I will also have tests for connectivity table

##### Surface Properties ######

@info "Tests for surface properties"

using JLD
data = load(joinpath(dirname(@__FILE__),"sphere.jld"))
#data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

n = Array{Float64}(undef,size(points)...)
NormalVectors!(n,points,faces,i->FaceVRing(i,faces))
curvatures = Array{Float64}(undef,size(points,2))
MeanCurvatures!(curvatures,points,faces,n,i->FaceVRing(i,faces))

angle_ = 0
curvaturer = 0

for xkey in 1:size(points,2)
    iter = FaceVRing(xkey,faces)
    ncalc = n[:,xkey]
    ccalc = curvatures[xkey]

    nT = points[:,xkey]/norm(points[:,xkey])
    angl = acos(dot(ncalc,nT))
    if angl>angle_
        global angle_ = angl
    end

    cT = 1.
    N = size(points,2)
    global curvaturer += abs((ccalc - cT)/cT)/N
end

@test angle_*180/pi < 1.3
@test curvaturer < 0.01

@test isapprox(volume(points,faces),4/3*pi,atol=0.05)
@test isapprox(*(FitEllipsoid(points)...),1)

### And now integrators
@info "Testing pasive mesh stabilisation"

function velocity(t,pos)
    x,y,z = pos
    
    x = x*0.15 + 0.35
    y = y*0.15 + 0.35
    z = z*0.15 + 0.35

    u = 2*sin(pi*x)^2 * sin(2*pi*y) * sin(2*pi*z) * sin(2/3*pi*t)
    v = - sin(2*pi*x) * sin(pi*y)^2 * sin(2*pi*z) * sin(2/3*pi*t)
    w = - sin(2*pi*x) * sin(2*pi*y) * sin(pi*z)^2 * sin(2/3*pi*t)

    [u,v,w] /0.15
end

v = zero(points)
t = 1
for i in 1:size(points,2)
    v[:,i] = velocity(t,points[:,i])
end

res =  stabilise(points,faces,v)
println("Energy before minimization Finit=$(res.Finit) after Fres=$(res.Fres)")



