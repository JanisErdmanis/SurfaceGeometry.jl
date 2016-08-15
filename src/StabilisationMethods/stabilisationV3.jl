using SurfaceGeometry
using Optim

using Parameters

@with_kw immutable Erdmanis2016
    C::AbstractFloat = 0.4
    ftol::AbstractFloat = 1e-6
end

function F(points,faces,v,zc::Erdmanis2016)

    Cp = zc.C
    
    s = 0    
    for ti in 1:size(faces,2)
        v1,v2,v3 = faces[:,ti]

        a = norm(points[:,v2] - points[:,v1])
        b = norm(points[:,v3] - points[:,v2])
        c = norm(points[:,v1] - points[:,v3])

        xava = dot(points[:,v2]-points[:,v1],v[:,v2]-v[:,v1])
        xbvb = dot(points[:,v3]-points[:,v2],v[:,v3]-v[:,v2])
        xcvc = dot(points[:,v1]-points[:,v3],v[:,v1]-v[:,v3])

        ### This part is responsible for first part in the sum
        s += v2>v1 ? 4 * xava^2 : 0
        s += v3>v2 ? 4 * xbvb^2 : 0
        s += v1>v3 ? 4 * xcvc^2 : 0

        ### This part is needed for ordinary things
        Cdelta = 1/4*sqrt(1 - 2* (a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        A = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( a^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        B = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( b^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        C = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( c^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        DCdelta = -A*xava - B*xbvb - C*xcvc

        s += Cp*DCdelta^2/Cdelta^2
    end

    return s
end

function gradF!(points,faces,v,storage,zc::Erdmanis2016)

    Cp = zc.C

    storage[:,:] = 0
    
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        a = norm(points[:,v2] - points[:,v1])
        b = norm(points[:,v3] - points[:,v2])
        c = norm(points[:,v1] - points[:,v3])

        Cdelta = 1/4*sqrt(1 - 2 * (a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        A = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( a^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        B = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( b^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        C = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( c^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        
        xa = points[:,v2] - points[:,v1]
        xb = points[:,v3] - points[:,v2]
        xc = points[:,v1] - points[:,v3]

        xava = dot(points[:,v2]-points[:,v1],v[:,v2]-v[:,v1])
        xbvb = dot(points[:,v3]-points[:,v2],v[:,v3]-v[:,v2])
        xcvc = dot(points[:,v1]-points[:,v3],v[:,v1]-v[:,v3])

        DCdelta = -A*xava - B*xbvb - C*xcvc

        A_ = 2*Cp*DCdelta*A/Cdelta^2 
        B_ = 2*Cp*DCdelta*B/Cdelta^2 
        C_ = 2*Cp*DCdelta*C/Cdelta^2 
        
        storage[:,v1] += (A_ - 8*xava)*xa - C_*xc
        storage[:,v2] += (B_ - 8*xbvb)*xb - A_*xa
        storage[:,v3] += (C_ - 8*xcvc)*xc - B_*xb
    end
end

function stabilise!(points,faces,n,v,zc::Erdmanis2016; op=:optim)
    stabiliseV3Optim!(points,faces,n,v,zc)
end

function stabiliseV3Optim!(points,faces,n,v,zc::Erdmanis2016)
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

    res = optimize(f,g!,v[:],method=Optim.ConjugateGradient(),ftol=zc.ftol)
    v[:,:] = reshape(res.minimum, size(v)...)[:,:]
end


function stabiliseV3!(points,faces,n,v,zc::Erdmanis2016)

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

