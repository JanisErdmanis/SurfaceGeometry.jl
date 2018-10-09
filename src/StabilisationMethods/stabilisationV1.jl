#using SurfaceGeometry

struct Zinchenko1997
    VertexRing
end

# function Zinchenko1997(points,faces,n;VertexRing=i -> VertexRing(i,faces))
#     Zinchenko1997(VertexRing)
# end    

eps = 1e-6

function f(i::Int,ni::Array{Float64,1},points,v,ring)
    s = zeros(3)
    pi = points[:,i]
    vi = v[:,i]
    for j in ring
        xij = points[:,j] - pi
        vij = v[:,j] - vi
        s += dot(xij,vij)*xij
    end
    P = eye(3) - ni*ni'
    return P*s
end

function fvec(i,points,v,ring)
    s = zeros(3)
    pi = points[:,i]
    vi = v[:,i]
    for j in ring
        xij = points[:,j] - pi
        vij = v[:,j] - vi
        s += dot(xij,vij)*xij
    end
    return s
end

function F(points,faces,v,mc::Zinchenko1997)

    function f(v1,v2)
        xij = points[:,v1] - points[:,v2]
        vij = v[:,v1] - v[:,v2]
        return dot(xij,vij)^2
    end
    
    N = size(faces,2)
    s = 0 
    for i in 1:N
        v1,v2,v3 = faces[:,i]
        s += v2>v1 ? f(v1,v2) : 0
        s += v3>v2 ? f(v2,v3) : 0
        s += v3>v1 ? f(v1,v3) : 0
    end
    return 4*s
end

# Returns ksi and eta for conjugate gradient iterations
function conjugrad(points,faces,va,vb,f)
    function A(i,j)
        xij = points[:,j] - points[:,i]
        vij = vb[:,j] - vb[:,i]
        return dot(xij,vij)
    end

    function B(i,j)
        xij = points[:,j] - points[:,i]
        fij = f[:,j] - f[:,i]
        return dot(xij,fij)
    end

    function C(i,j)
        xij = points[:,j] - points[:,i]
        Dvij = vb[:,j] - va[:,j] - (vb[:,i] - va[:,i])
        return dot(xij,Dvij)
    end

    S = zeros(2,2)
    P = zeros(2)

    function add!(i,j)
        a = A(i,j)
        b = B(i,j)
        c = C(i,j)

        S[1,1] += b^2
        S[1,2] += b*c
        S[2,1] += b*c
        S[2,2] += c^2

        P[1] += -a*b
        P[2] += -a*c
    end

    
    for i in 1:size(points,2)
        v1,v2,v3 = faces[:,i]

        if v2>v1
            add!(v1,v2)
        end

        if v3>v2
            add!(v3,v2)
        end

        if v3>v1
            add!(v1,v3)
        end
    end

    ksi, eta = S \ P
    return ksi,eta
end

# A simple conjugate gradient
function conjugrad0(points,faces,v,f)
    function A(i,j)
        xij = points[:,j] - points[:,i]
        vij = v[:,j] - v[:,i]
        return dot(xij,vij)
    end

    function B(i,j)
        xij = points[:,j] - points[:,i]
        fij = f[:,j] - f[:,i]
        return dot(xij,fij)
    end

    P = 0
    S = 0

    for i in 1:size(points,2)
        v1,v2,v3 = faces[:,i]

        if v2>v1
            a = A(v1,v2)
            b = B(v1,v2)
            S += b^2
            P += - a*b
        end

        
        if v3>v2
            a = A(v2,v3)
            b = B(v2,v3)
            S += b^2
            P += - a*b
        end

        if v3>v1
            a = A(v1,v3)
            b = B(v1,v3)
            S += b^2
            P += - a*b
        end
    end

    return P/S
end

function stabilise!(points,faces,n,va,mc::Zinchenko1997)

    VertexRing = mc.VertexRing
    
    ff = Array{Float64}(undef,size(points)...)
    for j in 1:size(points,2)
        ring = VertexRing(j)
        ff[:,j] = f(j,n[:,j],points,va,ring)
    end
    ksi0 = conjugrad0(points,faces,va,ff)
    vb = va + ksi0*ff
    vtemp = Array{Float64}(undef,size(va)...)

    Fitemp = Inf
    
    for i in 1:200
        
        for j in 1:size(points,2)
            ring = VertexRing(j)
            ff[:,j] = f(j,n[:,j],points,vb,ring)
        end

        ksi,eta = conjugrad(points,faces,va,vb,ff)

        vtemp[:,:] = vb[:,:]
        for j in 1:size(points,2)
            vb[:,j] = vb[:,j] + ksi*ff[:,j] + eta*(vb[:,j]-va[:,j])
        end
        va[:,:] = vtemp[:,:]

        Fi = F(points,faces,vb,mc)
        #println("F is $Fi and ksi is $ksi")

        if abs(Fitemp-Fi)<eps*Fi
            break

            if i>100
                warn("Stabilisation had $i iterations")
            end
        end
        Fitemp = Fi
    end
end
