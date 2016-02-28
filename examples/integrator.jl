using SurfaceGeometry

isdefined(:N) || (N=100)
isdefined(:PasiveStabilisation) || (PasiveStabilisation=false)
isdefined(:NormalVelocity) || (NormalVelocity=false)
isdefined(:ActiveStabilisation) || (ActiveStabilisation=false)
isdefined(:ViewerActive) || (ViewerActive=false)

using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]


function velocity(t,pos)
    x,y,z = pos
    
    x = x*0.15 + 0.35
    y = y*0.15 + 0.35
    z = z*0.15 + 0.35

    u = 2*sin(pi*x)^2 * sin(2*pi*y) * sin(2*pi*z) * sin(2/3*pi*t)
    v = - sin(2*pi*x) * sin(pi*y)^2 * sin(2*pi*z) * sin(2/3*pi*t)
    w = - sin(2*pi*x) * sin(2*pi*y) * sin(pi*z)^2 * sin(2/3*pi*t)

    [u,v,w] / 0.5 # / 0.15
end

### In reality this is very complicated function
function velocity{T<:Real}(t::Real,points::Array{T,2},faces::Array{Int,2})

    ### Some magic here on how to calculate the velocity
    v = Array(Float64,size(points)...)
    for xkey in 1:size(points,2)
        v[:,xkey] = velocity(t,points[:,xkey])
    end

    return v
end

memory = []

h = 3/N

t1 = 0
p1 = points
push!(memory,(t1,p1,faces))

### MidPoint method
v1 = velocity(t1,p1,faces)
v1 = velocity(t1 + h/2,p1 + h/2*v1,faces)
t2,p2 = t1 + h, p1 + h*v1

push!(memory,(t2,p2,faces))

for i in 1:50
    println("It is step $i")

    v2 = velocity(t2,p2,faces)
    
    ### Some stabilisation magic here
    # n = zeros(size(points)...)
    # NormalVectors!(n,points,faces,i->FaceVRing(i,faces))

    # for xkey in 1:size(points,2)
    #     #P = enye(3) - n[:,xkey]*n[:,xkey]'
    #     v2[:,xkey] = n[:,xkey]*dot(n[:,xkey],v2[:,xkey]) 
    # end

    if PasiveStabilisation==true
        res = stabilise(p2,faces,v2;vp=v1)
        println("Finit $(res.Finit) \t Fres $(res.Fres)")
        v2 = res.vres
    end

    ### Adams Bachforth oneliner
    (t1,p1,v1), (t2,p2) = (t2,p2,v2), (t2+h,p2 + h/2/(t2-t1)*( (2*(t2-t1) + h)*v2 - h*v1))

    push!(memory,(t2,p2,faces))
end


### For viewing results

using Escher
if !ViewerActive
    include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
    @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
    ViewerActive = true
end















