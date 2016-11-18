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
info("Mesh loading test")
using JLD
data = load(joinpath(dirname(@__FILE__),"sphere.jld"))
points = data["points"]
faces = data["faces"]

@test SphereError(points) < 0.01

info("Mesh generation with CGAL")

### Interface will change to pass signed distance function
### CGAL mesh generator

mesher = CGALSurfaceMesher()
fdis(x,y,z) = x^2 + y^2 + z^2 - 1
points, faces = SurfaceMesh(fdis,mesher)
@test SphereError(points) < 0.01

### Matlab mesh generator
info("Mesh generation with distmesh")
try
    step,boxsize = 0.2,[-1.1 -1.1 -1.1; 1.1 1.1 1.1]
    mesher = DistmeshSurfaceMesher(step,boxsize)

    a,b,c = 1,1,1
    signedf = """
    function f = fd(p)
        f = p(:,1).^2/$(a^2)+p(:,2).^2/$(b^2)+p(:,3).^2/$(c^2)-1;
    end
    
    """
    
    points, faces = SurfaceMesh(signedf,mesher)
    @test SphereError(points) < 0.01
catch
     warn("does not work")
end

### Testing a pushback
# sdist(x) = x[1]^2 + x[2]^2 + x[3]^2 - 1
# @test isapprox(norm(pushback(sdist,[1.,1.,1.])),1)

### Testing topology

info("Topology function tests")

t = [5,7,10,7,5,6,4,0,3,0,4,6,4,7,6,4,9,10,7,4,10,0,2,1,2,0,6,2,5,1,5,2,6,8,4,3,4,11,9,8,11,4,9,11,3,11,8,3]
faces = reshape(t,(3,div(length(t),3))) + 1 

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

info("Testing Complex DS")

points = zeros(size(faces)...)
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

info("Topology tests for Connectivity table")

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

info("Tests for surface properties")

using JLD
data = load(joinpath(dirname(@__FILE__),"sphere.jld"))
#data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

n = Array(Float64,size(points)...)
NormalVectors!(n,points,faces,i->FaceVRing(i,faces))
curvatures = Array(Float64,size(points,2))
MeanCurvatures!(curvatures,points,faces,n,i->FaceVRing(i,faces))

angle = 0
curvaturer = 0

for xkey in 1:size(points,2)
    iter = FaceVRing(xkey,faces)
    ncalc = n[:,xkey]
    ccalc = curvatures[xkey]

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
@test isapprox(*(FitEllipsoid(points)...),1)

### And now integrators

# Calculating velocities

mesher = CGALSurfaceMesher()
fdis(x,y,z) = x^2 + y^2 + z^2 - 1
points,faces = SurfaceMesh(fdis,mesher)

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

v = zeros(points)
t = 1
for i in 1:size(points,2)
    v[:,i] = velocity(t,points[:,i])
end

info("Testing pasive mesh stabilisation")
res =  stabilise(points,faces,v)
println("Energy before minimization Finit=$(res.Finit) after Fres=$(res.Fres)")

info("Testing ElTopo stabilisation")
par = Elparameters()
points, faces = improvemeshcol(points,faces,points + 0.1*v,par)

