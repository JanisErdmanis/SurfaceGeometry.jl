function halfedges!(vnormall,vects::AbstractArray{Float64,2})

    R = eye(3)
    #dump(vects)
    C,D,E = ZinchenkoDifertential!(vects,vnormall,R)

    z(x,y) = C*x^2 + D*x*y + E*y^2

    # halfvring = similar(vects)
    # println("From halfedges")
    Rinv = inv(R)
    # dump(Rinv)
    # dump(vects)
    for i in 1:size(vects,2)
        vec = vects[:,i] 
        vec = [vec[1]/2,vec[2]/2,z(vec[1]/2,vec[2]/2)]
        vec = Rinv*vec
        vects[:,i] = vec
    end
end

function smooth_surface_edge_nodes!(cmsh::TriGeom)
    
    points = cmsh.points
    faces = cmsh.faces
    enodes = cmsh.enodes
    epoints = cmsh.epoints
    epoints[:,:] = 0

    
    for v in 1:size(points,2)
        vring = Array(Int,0)
        for ti in TriRing(v,cmsh)
            face = cmsh.faces[:,ti]
            cw = (face.==v)[[3,1,2]]
            vi, = face[cw]
            push!(vring,vi)
        end

        p0 = points[:,v]
        vects = Array(Float64,3,length(vring))
        for (j,vi) in enumerate(vring)
            vects[:,j] = points[:,vi] - p0
        end

        ### Now updating positions

        vnormall = vnormal(v,cmsh)
        halfedges!(vnormall,vects)
        # ### points array must be set to zero for edge nodes
        
        ### Simple midpoint

        for (i,ti) in enumerate(TriRing(v,cmsh))
            face = cmsh.faces[:,ti]
            w = face .== v
            enode = cmsh.enodes[:,ti]
            enod, = enode[w[[2,3,1]]] ### Here might be an error

            value = vects[:,i] + p0 ### I can take average
            epoints[:,enod] += value/2
        end

    end

end    
