function find_other_triangle_edge{T<:Integer}(v1::Integer,v2::Integer,skip::Integer,t::AbstractArray{T,2})
    for i in 1:size(t)[end]
        if in(v1,t[:,i]) & in(v2,t[:,i]) & !(i==skip)
            return i
        end
    end
    return -1
end

function find_triangle_vertex{T<:Integer}(v::Integer,t::AbstractArray{T,2})

    for i in 1:size(t,2)
        if in(v,t[:,i])
            return i
        end
    end
    return size(t,2) # -1
end

immutable FaceRing
    v::Int
    faces::Array{Int,2}
end


function Base.start(iter::FaceRing)
    i = find_triangle_vertex(iter.v,iter.faces)
    return i
end

function Base.done(iter::FaceRing,i::Int)
    if i<size(iter.faces,2)
        return false
    else
        return true
    end
end

function Base.next(iter::FaceRing,i::Int)

    v = iter.v
    nexti = find_triangle_vertex(iter.v,iter.faces[:,i+1:end]) + i  # possible botleneck here
    return i, nexti
end


### As example I will use FBDS iterators already written

"""
Iterator for getting all triangles around vertex
"""
immutable TriRing
    v::Int
    faces::Array{Int,2}
    neighs::Array{Int,2} ### For now keeping simple
    vfaces::Array{Int,1}
    
    # function TriRing(v::Int,msh::TriGeom)

    #     faces = msh.faces
    #     neighs = msh.neighs
    #     vfaces = msh.vfaces
    #     new(v,faces,neighs,vfaces)
    # end

    # function TriRing(v::Int,faces::Array{Int,2},neighs::Array{Int,2},vfaces::Array{Int,1})
    #     new(v,faces,neighs,vfaces)
    # end
end



function Base.start(iter::TriRing)
    vface = iter.vfaces[iter.v]
    i0 = 1
    return (i0,vface)
end

function Base.done(iter::TriRing,state::Tuple{Int,Int})
    i, face = state
    i0, vface = start(iter)
    if !(i==i0) & (face==vface)
        return true
    else
        return false
    end    
end

function Base.next(iter::TriRing,state::Tuple{Int,Int})

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




