using SurfaceGeometry

isdefined(:N) || (N=100)
isdefined(:PasiveStabilisation) || (PasiveStabilisation=true)
isdefined(:TopologyStabilisation) || (TopologyStabilisation=true)
isdefined (:NormalVelocities) || (NormalVelocities=true)
isdefined(:ViewerActive) || (ViewerActive=false)

using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

function interpolate(p1,v1,p2)
    v2 = zeros(p2)
    for xkey in 1:size(p2,2)
        d = Inf
        vi = 0
        for ykey in 1:size(p1,2)
            dist = norm(p1[:,ykey] - p2[:,xkey])
            if dist<d
                d = dist
                vi = ykey
            end
        end
        v2[:,xkey] = v1[:,vi]
    end
end

function velocity(t,points,faces) 
    v = Array(Float64,size(points)...)
    if NormalVelocities==true
        normals = Array(Float64,size(points)...);
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    

        for xkey in 1:size(points,2)
            v[:,xkey] = velocity(t,points[:,xkey])
            nx = normals[:,xkey]
            v[:,xkey] = nx*dot(nx,v[:,xkey])
        end
    else
        for xkey in 1:size(points,2)
            v[:,xkey] = velocity(t,points[:,xkey])
        end
    end
    return v
end

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

memory = []

h = 3/N

t = 0.001
p = points
f = faces
push!(memory,(t,p,f))

vp = zeros(p)
pp = zeros(p)

v1 = zeros(p)
v2 = zeros(p)

zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)

scale = 0.2
par = Elparameters(
                   m_use_fraction = false,
                   m_min_edge_length = 0.7*scale,
                   m_max_edge_length = 1.5*scale,
                   m_max_volume_change = 0.1*scale^3,
                   m_min_curvature_multiplier = 1.0,
                   m_max_curvature_multiplier = 1.0,
                   m_merge_proximity_epsilon = 0.5*scale,
                   m_proximity_epsilon = 0.00001,
m_perform_improvement = true, 
m_collision_safety = false,
m_min_triangle_angle = 15,
m_max_triangle_angle = 120,
m_allow_vertex_movement = true, ### 
m_use_curvature_when_collapsing = false,
m_use_curvature_when_splitting = false,
m_dt = h
)

for i in 1:35
    println("step $i just started")
    
    v1 = velocity(t,p,f)
    if PasiveStabilisation==true
        res = stabilise(p,f,v1;vp=interpolate(pp,vp,p),method=zc)
        println("Finit $(res.Finit) \t Fres $(res.Fres)")
        v1 = res.vres
    end

    v2 = velocity(t + h/2,p + h/2*v1,f)
    if PasiveStabilisation==true
        res = stabilise(p + h/2*v1,f,v2;vp=v1,method=zc)
        println("Finit $(res.Finit) \t Fres $(res.Fres)")
        v2 = res.vres
    end

    vp = copy(v2)
    pp = copy(p)
    
    ### Adams Bachforth oneliner
    #(t1,p1,v1), (t2,p2) = (t2,p2,v2), (t2+h,p2 + h/2/(t2-t1)*( (2*(t2-t1) + h)*v2 - h*v1))

    if TopologyStabilisation==true
        actualdt,p,f = improvemeshcol(p,f,p + h*v2,par)
        t = t + actualdt
    else
        p = p + h*v2
        t = t + h
    end        

    push!(memory,(t,p,f))
end

using Escher
if !ViewerActive
    include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
    @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
    ViewerActive = true
end
