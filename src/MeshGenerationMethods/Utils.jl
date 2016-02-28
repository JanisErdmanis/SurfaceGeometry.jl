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

