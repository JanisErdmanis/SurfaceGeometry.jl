### More advanced data structures

struct FaceBasedDS
    faces::AbstractArray{Int,2}
    neighs::AbstractArray{Int,2} ### For now keeping simple
    vfaces::AbstractArray{Int,1}
end


function FaceBasedDS(faces::Array{Int,2})

    vfaces = Array{Int}(undef,maximum(faces))
    neighs = Array{Int}(undef,size(faces)...)
    
    for vi in 1:maximum(faces)
        vfaces[vi] = find_triangle_vertex(vi,faces)
    end

    for ti in 1:size(faces,2)
        v1, v2, v3 = faces[:,ti]
        t1 = find_other_triangle_edge(v2,v3,ti,faces)
        t2 = find_other_triangle_edge(v1,v3,ti,faces)
        t3 = find_other_triangle_edge(v1,v2,ti,faces)
        neighs[:,ti] = [t1,t2,t3]
    end
    
    return FaceBasedDS(faces,neighs,vfaces)
end

"""
Iterator for getting all triangles around vertex
"""
struct FaceVRingFB
    v::Int
    faces::Array{Int,2}
    neighs::Array{Int,2} ### For now keeping simple
    vfaces::Array{Int,1}
end

FaceVRing(v::Int,ds::FaceBasedDS) = FaceVRingFB(v,ds.faces,ds.neighs,ds.vfaces)
    
function start(iter::FaceVRingFB)
    vface = iter.vfaces[iter.v]
    i0 = 1
    return (i0,vface)
end

function done(iter::FaceVRingFB,state::Tuple{Int,Int})
    i, face = state
    i0, vface = start(iter)
    if !(i==i0) & (face==vface)
        return true
    else
        return false
    end    
end

function next(iter::FaceVRingFB,state::Tuple{Int,Int})

    i, tri = state
    v = iter.v
    face = iter.faces[:,tri]
    neighbours = iter.neighs[:,tri]

    index = face.==v
    w = index[[1,2,3]]
    cw = index[[3,1,2]]

    nexttri, = neighbours[cw]

    if nexttri==-1
        error("The surface is not closed")
    end
    
    return tri, (i+1,nexttri)
end


Base.iterate(iter::FaceVRingFB) = next(iter,start(iter))
function Base.iterate(iter::FaceVRingFB,ti)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end


###

struct DoubleVertexVRingFB
    v::Int
    faces::Array{Int,2}
    neighs::Array{Int,2} ### For now keeping simple
    vfaces::Array{Int,1}
end

DoubleVertexVRing(v::Int,ds::FaceBasedDS) = DoubleVertexVRingFB(v,ds.faces,ds.neighs,ds.vfaces)

function start(iter::DoubleVertexVRingFB)
    vface = iter.vfaces[iter.v]
    i0 = 1
    return (i0,vface)
end

function done(iter::DoubleVertexVRingFB,state::Tuple{Int,Int})
    i, face = state
    i0, vface = start(iter)
    if !(i==i0) & (face==vface)
        return true
    else
        return false
    end    
end

function next(iter::DoubleVertexVRingFB,state::Tuple{Int,Int})

    i, tri = state
    v = iter.v
    face = iter.faces[:,tri]
    neighbours = iter.neighs[:,tri]

    index = face.==v
    w = index[[1,2,3]]
    cw = index[[3,1,2]]
    ccw = index[[2,3,1]]

    nexttri, = neighbours[cw]

    if nexttri==-1
        error("The surface is not closed")
    end
    
    return (face[cw]...,face[ccw]...), (i+1,nexttri)
end

Base.iterate(iter::DoubleVertexVRingFB) = next(iter,start(iter))
function Base.iterate(iter::DoubleVertexVRingFB,ti)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end


"""
Now the magic of making VertexIterator
"""
struct VertexVRingFB
    v::Int
    faces::Array{Int,2}
    neighs::Array{Int,2} ### For now keeping simple
    vface::Int
    # function VertexRingFBDS(v,faces,neighs,vfaces)
    #     new(v,faces,neighs,vfaces[v])
    # end
end

VertexVRing(v::Int,ds::FaceBasedDS) = VertexVRingFB(v,ds.faces,ds.neighs,ds.vfaces[v])

function start(iter::VertexVRingFB)
    vface = iter.vface
    i0 = 1
    return (i0,vface)
end

function done(iter::VertexVRingFB,state::Tuple{Int,Int})
    i, face = state
    if !(i==1) & (face==iter.vface)
        return true
    else
        return false
    end    
end

function next(iter::VertexVRingFB,state::Tuple{Int,Int})

    i, tri = state
    v = iter.v
    face = iter.faces[:,tri]
    neighbours = iter.neighs[:,tri]

    index = face.==v
    w = index[[1,2,3]]
    cw = index[[3,1,2]]

    nexttri, = neighbours[cw]

    if nexttri==-1
        error("The surface is not closed")
    end

    ### Code for extracting vertex from face tri

    face = iter.faces[:,tri]
    cw = (face.==v)[[3,1,2]]
    vi, = face[cw]
    
    return vi, (i+1,nexttri)
end

Base.iterate(iter::VertexVRingFB) = next(iter,start(iter))
function Base.iterate(iter::VertexVRingFB,ti)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end

### Using conectivity for improving performance

struct ConnectivityDS
    connectivity
end

function ConnectivityTable(faces,valence)
    vmax = maximum(faces)
    connectivity = zeros(Int,valence,vmax)
    for i in 1:vmax
        v1 = Array{Int}(undef,0)
        v2 = Array{Int}(undef,0)
        for (v1i,v2i) in DoubleVertexVRing(i,faces)
            push!(v1,v1i)
            push!(v2,v2i)
        end
        
        vj = v2[1]
        connectivity[1,i] = vj
        for j in 2:size(v1,1)
            w = (v1.==vj)
            vj, = v2[w]
            connectivity[j,i] = vj
        end
    end
    return connectivity
end

ConnectivityDS(faces,valence) = ConnectivityDS(ConnectivityTable(faces,valence))

struct VertexVRingCon{T<:Integer}
    v::T
    connectivity::Array{T,2}
end

VertexVRing(v::Int,ds::ConnectivityDS) = VertexVRingCon(v,ds.connectivity)

function start(iter::VertexVRingCon)
    return 1
end

function done(iter::VertexVRingCon,i::T) where {T<:Integer}
    return (i>size(iter.connectivity,1)) || (iter.connectivity[i,iter.v]==0) ? true : false
end

function next(iter::VertexVRingCon,i::T) where {T<:Integer}
    return iter.connectivity[i,iter.v], i+1
end


Base.iterate(iter::VertexVRingCon) = next(iter,start(iter))
function Base.iterate(iter::VertexVRingCon,ti)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end


struct DoubleVertexVRingCon{T<:Integer}
    v::T
    connectivity::Array{T,2}
end

DoubleVertexVRing(v::Int,ds::ConnectivityDS) = DoubleVertexVRingCon(v,ds.connectivity)

function start(iter::DoubleVertexVRingCon)
    return 1
end

function done(iter::DoubleVertexVRingCon,i::T) where {T<:Integer}
    return (i>size(iter.connectivity,1)) || (iter.connectivity[i,iter.v]==0) ? true : false
end

function next(iter::DoubleVertexVRingCon,i::T) where {T<:Integer}
    v1 = iter.connectivity[i,iter.v]
    v2 = done(iter,i+1) ? iter.connectivity[1,iter.v] : iter.connectivity[i+1,iter.v] # (i==size(iter.connectivity,1)) || (iter.connectivity[i+1,iter.v]==0)
    return (v1,v2), i+1
end

Base.iterate(iter::DoubleVertexVRingCon) = next(iter,start(iter))
function Base.iterate(iter::DoubleVertexVRingCon,ti)
    if done(iter,ti)
        return nothing
    else
        return next(iter,ti)
    end
end
