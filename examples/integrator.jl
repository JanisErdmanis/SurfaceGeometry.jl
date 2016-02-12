using SurfaceGeometry

isdefined(:N) || (N=100)
isdefined(:meth) || (meth=AdamsStep)
isdefined(:PasiveStabilisation) || (PasiveStabilisation=false)
isdefined(:ActiveStabilisation) || (ActiveStabilisation=false)
isdefined(:ViewerActive) || (ViewerActive=false)

using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

function velocity(points::Array{Float64,2},faces::Array{Int,2},t::Float64) ### Iterator also can be fit as an argument
    v = Array(Float64,size(points)...)
    for xkey in 1:size(points,2)
        v[:,xkey] = velocity(points[:,xkey],t)
    end
    return v
end

function velocity(pos,t)
    x,y,z = pos
    
    x = x*0.15 + 0.35
    y = y*0.15 + 0.35
    z = z*0.15 + 0.35

    u = 2*sin(pi*x)^2 * sin(2*pi*y) * sin(2*pi*z) * sin(2/3*pi*t)
    v = - sin(2*pi*x) * sin(pi*y)^2 * sin(2*pi*z) * sin(2/3*pi*t)
    w = - sin(2*pi*x) * sin(2*pi*y) * sin(pi*z)^2 * sin(2/3*pi*t)

    [u,v,w]  # /0.15
end

steppers = convert(Array{meth,1},[])
for i in 1:size(points,2)
    s = meth(0,points[:,i])
    push!(steppers,s)
end

ti = 0.
points = step!(steppers,velocity(points,faces,ti),0.000001)
ti += 0.000001

h = 3/N

t = []
frames = []
push!(t,ti)
push!(t,(points[:,:],faces[:,:]))    

for i in 1:N
    v = velocity(points,faces,ti)

    if PasiveStabilisation==true
        for xkey in 1:size(points,2)
            iter = FaceRing(xkey,faces)
            n = vnormal(xkey,points,faces,iter)
            v[:,xkey] = n*dot(n,v[:,xkey])
        end
    end
    
    points = step!(steppers,v,h)
    ti += h
    push!(t,ti)
    push!(frames,(points[:,:],faces[:,:]))    
end

### For viewing results

using Escher
if !ViewerActive
    include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
    @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
    ViewerActive = true
end

