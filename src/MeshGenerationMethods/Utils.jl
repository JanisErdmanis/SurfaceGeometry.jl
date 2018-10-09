# Shared methods between surface meshers
# * pushback
# * isoriented
# * volume (because isoriented could be made more robust)
# * subdivision


"Calculation of precise volume. Needed for convergence check"
function volume(points,faces)

    normal0 = [0,0,1]

    s = 0

    for tri in 1:size(faces,2)
        face = faces[:,tri]
        y1 = points[:,face[1]]
        y2 = points[:,face[2]]
        y3 = points[:,face[3]]

        normaly = cross(y2-y1,y3-y1)
        normaly /= norm(normaly)

        area = norm(cross(y2-y1,y3-y1))/2
        areaproj = dot(normaly,normal0)*area
        volume = dot(y1 + y2 + y3,normal0)/3*areaproj

        s += volume
    end

    return s

    # Calculate face normal
    # Calculate projected area
    # Calculate ordinary volume 
    # Calculate volume between projected and real area
    # (+) if normal is outwards
end


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


function isoriented(points::Array{T,2},triangles::Array{S,2}) where {T <: AbstractFloat, S <: Integer}

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


function subdivision(faces)
    edges = []

    for ti in 1:size(faces,2)
        v1,v2,v3 = faces[:,ti]

        v2<v1 || push!(edges,(v1,v2))
        v3<v2 || push!(edges,(v2,v3))
        v1<v3 || push!(edges,(v3,v1))
    end


    rfaces = Array(Int,3,size(faces,2)*4)
    N = maximum(faces)

    for ti in 1:size(faces,2)

        v1,v2,v3 = faces[:,ti]
        e3 = findfirst(edges,v2>v1 ? (v1,v2) : (v2,v1)) + N 
        e1 = findfirst(edges,v3>v2 ? (v2,v3) : (v3,v2)) + N 
        e2 = findfirst(edges,v1>v3 ? (v3,v1) : (v1,v3)) + N 

        rfaces[:,4*ti - 3] = [v1,e3,e2]
        rfaces[:,4*ti - 2] = [v2,e1,e3]
        rfaces[:,4*ti - 1] = [v3,e2,e1]
        rfaces[:,4*ti] = [e1,e2,e3]
    end

    return rfaces    
end

function subdivision(points,faces; n=1,method=:linear)
    if method==:linear
        rpoints,rfaces = LinearSubdivision(points,faces)
    elseif method==:paraboloid
        rpoints,rfaces = ParaboloidSubdivision(points,faces)
    end

    return rpoints,rfaces
end

function LinearSubdivision(points,faces)

    rfaces = subdivision(faces)
    rpoints = Array(Float64,3,maximum(rfaces))

    for ti in 1:size(faces,2)

        v1,v2,v3 = faces[:,ti]
        e1,e2,e3 = rfaces[:,4*ti]

        rpoints[:,v1] = points[:,v1]
        rpoints[:,v2] = points[:,v2]
        rpoints[:,v3] = points[:,v3]
        rpoints[:,e1] = (points[:,v2] + points[:,v3])/2
        rpoints[:,e2] = (points[:,v1] + points[:,v3])/2
        rpoints[:,e3] = (points[:,v1] + points[:,v2])/2
    end

    return rpoints,rfaces
end


function ParaboloidSubdivision(points,faces)

    normals = Array(Float64,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    
    rfaces = subdivision(faces)
    rpoints = Array(Float64,3,maximum(rfaces))

    vproperties = []
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]
        
        vects = Float64[]
        for i in VertexVRing(xkey,faces)
            vi = points[:,i] - x
            vects = [vects..., vi...]
        end
        vects = reshape(vects,3,div(length(vects),3))

        R = eye(3)
        C,D,E = SurfaceGeometry.ZinchenkoDifertential!(vects,nx,R)

        push!(vproperties,(C,D,E,R))
    end

    z(x,y,C,D,E) = C*x^2 + D*x*y + E*y^2
    pushon(x,C,D,E,R) = ((xp,yp,zp) = R*x; inv(R)*[xp,yp,z(xp,yp,C,D,E)])
    
    for ti in 1:size(faces,2)
        v1,v2,v3 = faces[:,ti]
        e1,e2,e3 = rfaces[:,4*ti]

        rpoints[:,v1] = points[:,v1]
        rpoints[:,v2] = points[:,v2]
        rpoints[:,v3] = points[:,v3]

        ### It gets overwritten once        
        xe1 = (points[:,v2] + points[:,v3])/2
        rpoints[:,e1] = xe1 + (pushon(xe1-points[:,v2],vproperties[v2]...) + pushon(xe1-points[:,v3],vproperties[v3]...))/2

        xe2 = (points[:,v3] + points[:,v1])/2
        rpoints[:,e2] = xe2 + (pushon(xe2-points[:,v1],vproperties[v1]...) + pushon(xe2-points[:,v3],vproperties[v3]...))/2

        xe3 = (points[:,v1] + points[:,v2])/2
        rpoints[:,e3] = xe3 + (pushon(xe3-points[:,v1],vproperties[v1]...) + pushon(xe3-points[:,v2],vproperties[v2]...))/2
    end

    return rpoints,rfaces
end

"Subdivision of surface and pushback"
function subdivision(points,faces,sdist::Function)
    rpoints,rfaces = LinearSubdivision(points,faces)

    for i in 1:size(rpoints,2)
        rpoints[:,i] = pushback(sdist,rpoints[:,i])
    end

    return rpoints,rfaces
end
