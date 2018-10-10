function find_other_triangle_edge(v1::Integer,v2::Integer,skip::Integer,t::AbstractArray{T,2}) where {T<:Integer}
    for i in 1:size(t)[end]
        if in(v1,t[:,i]) & in(v2,t[:,i]) & !(i==skip)
            return i
        end
    end
    return -1
end

function find_triangle_vertex(v::Integer,t::AbstractArray{T,2}) where {T<:Integer}

    for i in 1:size(t,2)
        if in(v,t[:,i])
            return i
        end
    end
    return size(t,2) + 1 # -1
end

### Deprecaction Warnings
FaceRing(v,faces) = error("Deprecated: use FaceVRing")
FaceRingS(v,faces) = error("Deprecated: use DoubleVertexVRing")
VeretexRing(v,faces) = error("Deprecated: use VertexVRing")

struct FaceVRing
    v::Int
    faces::Array{Int,2}
end

start(iter::FaceVRing) = find_triangle_vertex(iter.v,iter.faces)
done(iter::FaceVRing,ti::Int) = ti<=size(iter.faces,2) ? false : true

function next(iter::FaceVRing,i::Int)

     v = iter.v
    nexti = find_triangle_vertex(iter.v,iter.faces[:,i+1:end]) + i  # possible botleneck here
    return i, nexti
end

Base.iterate(iter::FaceVRing) = next(iter,start(iter))
function Base.iterate(iter::FaceVRing,ti::Int)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end

struct DoubleVertexVRing
    v::Int
    faces::Array{Int,2}
end

start(iter::DoubleVertexVRing) = find_triangle_vertex(iter.v,iter.faces)
done(iter::DoubleVertexVRing,ti::Int) = ti<=size(iter.faces,2) ? false : true

function next(iter::DoubleVertexVRing,i::Int)

    v = iter.v
    nexti = find_triangle_vertex(iter.v,iter.faces[:,i+1:end]) + i  # possible botleneck here
    face = iter.faces[:,i]
    w = face.==v
    cw = w[[3,1,2]]
    ccw = w[[2,3,1]]
    
    return (face[cw]...,face[ccw]...), nexti
end

Base.iterate(iter::DoubleVertexVRing) = next(iter,start(iter))
function Base.iterate(iter::DoubleVertexVRing,ti::Int)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end

### Unordered Veretex ring

struct VertexVRing
    v::Int
    faces::Array{Int,2}
end

start(iter::VertexVRing) = find_triangle_vertex(iter.v,iter.faces)
done(iter::VertexVRing,ti::Int) = ti<=size(iter.faces,2) ? false : true

function next(iter::VertexVRing,ti::Int)

    v = iter.v
    faces = iter.faces
    face = faces[:,ti]
    cw = (face.==v)[[3,1,2]]
    vi, = face[cw]
        
    nexti = find_triangle_vertex(iter.v,iter.faces[:,ti+1:end]) + ti  # possible botleneck
    return vi, nexti
end

Base.iterate(iter::VertexVRing) = next(iter,start(iter))
function Base.iterate(iter::VertexVRing,ti::Int)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end












