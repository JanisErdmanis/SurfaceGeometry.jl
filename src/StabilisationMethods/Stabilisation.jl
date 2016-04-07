### All possible fields
### If it initialisation lacks some fields then it might be resolved with simpler constructors
### Similar type could be made for RestructurisationResults
immutable StabilisationResults
    zc
    vinit
    Finit
    vres
    Fres
end

function Base.show(io::IO, r::StabilisationResults)
    println("Results of Stabilisation Algorithm $(typeof(r.zc))")
    println(" * vinit: $(r.vinit)")
    println(" * Finit: $(r.Finit)")
    println(" * vres: $(r.vres)")
    println(" * Fres: $(r.Fres)")
end

# High level interface. If you trust enough ;)
function stabilise(points,faces,v;method=:Zinchenko2013,vp=zeros(size(points)...))

    n = Array(Float64,size(points)...)
    NormalVectors!(n,points,faces,i -> FaceVRing(i,faces))
    
    ### Projecting previous velocity
    for xkey in 1:size(points,2)
        P = eye(3) - n[:,xkey]*n[:,xkey]'
        v[:,xkey] = n[:,xkey]*dot(n[:,xkey],v[:,xkey]) + P*vp[:,xkey]
    end
    
    if method==:Zinchenko2013
        zc = Zinchenko2013(points,faces,n)
    elseif method==:Zinchenko1997
        zc = Zinchenko1997(i -> VertexVRing(i,faces))
    else
        zc = method   ### For construction outside the usual 
    end

    vres = copy(v)
    Finit = F(points,faces,v,zc)
    stabilise!(points,faces,n,vres,zc)
    Fres = F(points,faces,vres,zc)
    res = StabilisationResults(zc,v,Finit,vres,Fres)
    return res
end


# function stabilise(points,faces,v;method=:Zinchenko2013)
    
#     n = Array(Float64,size(points)...)
#     for i in 1:size(points,2)
#         iter = FaceVRing(i,faces)
#         n[:,i] = vnormal(i,points,faces,iter)
#     end

#     stabilise(points,faces,n,v;method=method)
# end



