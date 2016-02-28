using SurfaceGeometry
using Optim

immutable Zinchenko2013{T<:AbstractFloat}
    gamma::AbstractFloat 
    C::AbstractFloat
    h2::Vector{T}
end

function Zinchenko2013(points,faces,n;gamma=0.25,C=0.4)

    function VArea(vi)
        area = 0
        for tri in FaceVRing(vi,faces)
            v0,v1,v2 = faces[:,tri]
            y0 = points[:,v0]
            y1 = points[:,v1]
            y2 = points[:,v2]

            Svec = cross(y1-y0,y2-y0)/2
            triarea = norm(Svec)

            area += triarea
        end
        return area
    end

    Lambda = Array(Float64,size(points,2))
    for i in 1:size(points,2)

        normal = n[:,i]

        p0 = points[:,i]
        vects = Array(Float64,3,0)
        for vi in VertexVRing(i,faces)
            vects = [vects (points[:,vi]-p0)]
        end

        R = eye(3)
        C,D,E = ZinchenkoDifertential!(vects,normal,R)

        A = [C D/2;D/2 E]
        k1,k2 = eigvals(-A)*2

        lambda = k1^2 + k2^2
        Lambda[i] = lambda
    end

    s = 0
    #gamma = 0.25
    N = size(points,2)
    for i in 1:size(points,2)
        s += Lambda[i]^gamma * VArea(i)/3
    end
    K = 4/sqrt(3)/N * s

    h2 = Array(Float64,size(points,2))
    for i in 1:size(points,2)
        h2[i] = K*Lambda[i]^(-gamma)
    end

    Zinchenko2013(gamma,C,h2)
end


function F(points,faces,v,zc::Zinchenko2013)

    h2 = zc.h2
    
    s = 0    
    for ti in 1:size(faces,2)
        v1,v2,v3 = faces[:,ti]

        a = norm(points[:,v2] - points[:,v1])
        b = norm(points[:,v3] - points[:,v2])
        c = norm(points[:,v1] - points[:,v3])

        xava = dot(points[:,v2]-points[:,v1],v[:,v2]-v[:,v1])
        xbvb = dot(points[:,v3]-points[:,v2],v[:,v3]-v[:,v2])
        xcvc = dot(points[:,v1]-points[:,v3],v[:,v1]-v[:,v3])

        s += v2>v1 ? 4*(2/(h2[v2]+h2[v1]) - (h2[v2] + h2[v1])/2/a^4)^2 * xava^2 : 0
        s += v3>v2 ? 4*(2/(h2[v2]+h2[v3]) - (h2[v2] + h2[v3])/2/b^4)^2 * xbvb^2 : 0
        s += v1>v3 ? 4*(2/(h2[v1]+h2[v3]) - (h2[v1] + h2[v3])/2/c^4)^2 * xcvc^2 : 0
        
        Cdelta = 1/4*sqrt(1 - 2* (a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        A = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( a^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        B = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( b^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        C = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( c^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        DCdelta = -A*xava - B*xbvb - C*xcvc

        s += 0.4*DCdelta^2/Cdelta^2
    end

    return s
end

function gradF!(points,faces,v,storage,zc::Zinchenko2013)

    h2 = zc.h2

    storage[:,:] = 0
    
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        a = norm(points[:,v2] - points[:,v1])
        b = norm(points[:,v3] - points[:,v2])
        c = norm(points[:,v1] - points[:,v3])

        xava = dot(points[:,v2]-points[:,v1],v[:,v2]-v[:,v1])
        xbvb = dot(points[:,v3]-points[:,v2],v[:,v3]-v[:,v2])
        xcvc = dot(points[:,v1]-points[:,v3],v[:,v1]-v[:,v3])
        
        Cdelta = 1/4*sqrt(1 - 2 * (a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        A = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( a^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        B = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( b^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        C = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( c^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        DCdelta = -A*xava - B*xbvb - C*xcvc

        
        ha2 = (h2[v2] + h2[v1])/2
        hb2 = (h2[v3] + h2[v2])/2
        hc2 = (h2[v1] + h2[v3])/2
        
        A_ = 2*0.4*DCdelta*A/Cdelta^2 
        B_ = 2*0.4*DCdelta*B/Cdelta^2 
        C_ = 2*0.4*DCdelta*C/Cdelta^2 

        xa = points[:,v2] - points[:,v1]
        xb = points[:,v3] - points[:,v2]
        xc = points[:,v1] - points[:,v3]
        
        storage[:,v1] += (A_ - 8*(1/ha2 - ha2/a^4)^2*xava)*xa - C_*xc
        storage[:,v2] += (B_ - 8*(1/hb2 - hb2/b^4)^2*xbvb)*xb - A_*xa
        storage[:,v3] += (C_ - 8*(1/hc2 - hc2/c^4)^2*xcvc)*xc - B_*xb

    end
end

function stabilise!(points,faces,n,v,zc::Zinchenko2013; op=:optim)
    if op==:optim
        stabiliseV2Optim!(points,faces,n,v,zc)
    elseif op===:zinchenko
        stabiliseV2!(points,faces,n,v,zc)
    else
        error("Your requested method is not available")
    end
end

    

function stabiliseV2Optim!(points,faces,n,v,zc::Zinchenko2013)
    #h2 = hvec(points,faces,n)

    function f(x::Vector)
        vv = reshape(x,3,div(length(x),3))
        return F(points,faces,vv,zc)
    end


    function g!(x::Vector, storage::Vector)

        vv = reshape(x,size(points)...)
        gradf = reshape(storage,size(points)...)

        gradF!(points,faces,vv,gradf,zc)
        
        for i in 1:size(v,2)
            P = eye(3)-n[:,i]*n[:,i]'
            gradf[:,i] = P*gradf[:,i]
        end
    end

    res = optimize(f,g!,v[:],method=:cg,ftol=1e-6)
    v[:,:] = reshape(res.minimum, size(v)...)[:,:]
end


function stabiliseV2!(points,faces,n,v,zc::Zinchenko2013)

    gradF_v = Array(Float64,size(points)...)
    gradF!(points,faces,v,gradF_v,zc)
    
    f = Array(Float64,size(points)...)
    for i in 1:size(points,2)
        f[:,i] = - (eye(3) - n[:,i]*n[:,i]') * gradF_v[:,i] 
    end

    Fs = Inf
    eps = 1e-4

    gradF_f = Array(Float64,size(points)...)
    gradF!(points,faces,f,gradF_f,zc)

    q = 2*dot(gradF_f[:],f[:])
    s = dot(gradF_v[:],f[:]) + dot(gradF_f[:],v[:])
    ksi = - s/q

    gradF_dv = Array(Float64,size(points)...)
    vp = copy(v)
    v[:,:] = v[:,:] + ksi*f[:,:]

    step = 1
    for step in 1:150

        gradF!(points,faces,v,gradF_v,zc)
        for i in 1:size(points,2)
            f[:,i] = -(eye(3) - n[:,i]*n[:,i]') * gradF_v[:,i]  
        end

        gradF!(points,faces,f,gradF_f,zc)
        gradF!(points,faces,v-vp,gradF_dv,zc)
                
        S = [dot(gradF_v[:],f[:]) + dot(gradF_f[:],v[:]), dot(gradF_dv[:],v[:]) + dot(gradF_v[:],v[:] - vp[:])]
        Q = [2*dot(gradF_f[:],f[:])   dot(gradF_f[:],v[:]-vp[:])+dot(gradF_dv[:],f[:]); dot(gradF_dv[:],f[:])+dot(gradF_v[:],v[:]-vp[:])   2*dot(gradF_dv[:],v[:]-vp[:])]
        ksi, eta = -Q\S

        v[:,:],vp[:,:] = v[:,:] + ksi*f[:,:]  + eta*(v[:,:] - vp[:,:]), v[:,:]

        Fsnew = F(points,faces,v,zc)
        if (abs(Fs - Fsnew) < eps*Fsnew )
            break
        else
            Fs = Fsnew
        end
    end

    if step>100
        warn("stabilistation used $step steps")
    end
end

