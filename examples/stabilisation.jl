using SurfaceGeometry

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

ti = 1.
v = velocity(points,faces,ti)

FBDS = FaceBasedDS(faces)
n = Array(Float64,size(points)...)
for i in 1:size(points,2)
    #iter = FaceVRing(i,faces)
    iter = FaceVRing(i,FBDS)    
    n[:,i] = vnormal(i,points,faces,iter)
    v[:,i] = dot(n[:,i],v[:,i])*n[:,i]
end

#res = stabilise(points,faces,n,v)
res = stabilise(points,faces,n,v;method=:Zinchenko2013)



