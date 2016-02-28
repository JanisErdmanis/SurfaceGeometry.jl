function vnormal(v::Int,points::Array{Float64,2},faces::Array{Int,2},Iterator)

    N_MWA = zeros(3)

    for ti in Iterator
        face = faces[:,ti]

        w = face.==v
        cw = w[[3,1,2]]
        ccw = w[[2,3,1]]

        p1 = points[:,face[w]...]
        p2 = points[:,face[cw]...]
        p3 = points[:,face[ccw]...]

        a1 = p2 - p1
        a2 = p3 - p1

        alpha = asin(norm(cross(a1,a2))/norm(a1)/norm(a2))

        N_MWA += alpha*cross(a1,a2)
        
    end

    N_MWA = N_MWA/norm(N_MWA)

    return N_MWA
end

function NormalVectors!(n,points,faces,FaceVRing)
    for i in 1:size(points,2)
        n[:,i] = vnormal(i,points,faces,FaceVRing(i))
    end
end

function ZinchenkoDifertential!(vects,vnormal,R)

    Lx = [0 0 0; 0 0 -1; 0 1 0]
    Ly = [0 0 1; 0 0 0; -1 0 0]
    Lz = [0 -1 0; 1 0 0; 0 0 0]

    d = [0,0,1] + vnormal
    d /= norm(d)
        
    Ln = d[1]*Lx + d[2]*Ly + d[3]*Lz
    R[:,:] = expm(pi*Ln)

    # println("Rotation matrix is")
    # dump(R)

    for vj in 1:size(vects)[end]
        vects[:,vj] = R*vects[:,vj]
    end

    ### Construction of the system
    A = Array(Float64,3,3)
    B = Array(Float64,3)

    vects_norm2 = Array(Float64,size(vects)[end])
    for vj in 1:size(vects)[end]
       vects_norm2[vj] = norm(vects[:,vj])^2
    end

    A[1,1] = sum(vects[1,:].^4./vects_norm2)
    A[1,2] = sum(vects[1,:].^3.*vects[2,:]./vects_norm2)
    A[1,3] = sum(vects[2,:].^2.*vects[1,:].^2./vects_norm2)
    A[2,1] = A[1,2]
    A[2,2] = A[1,3]
    A[2,3] = sum(vects[2,:].^3.*vects[1,:]./vects_norm2)
    A[3,1] = A[1,3]
    A[3,2] = A[2,3]
    A[3,3] = sum(vects[2,:].^4./vects_norm2)

    B[1] = sum(vects[3,:].*vects[1,:].^2./vects_norm2)
    B[2] = sum(vects[1,:].*vects[2,:].*vects[3,:]./vects_norm2)
    B[3] = sum(vects[2,:].^2.*vects[3,:]./vects_norm2)

    C,D,E = A\B
    return C,D,E
end

### Orginals ieks utils.jl
function vcurvature(v::Int,points::Array{Float64,2},faces::Array{Int,2},normal,Iterator)

    #points = msh.points
    vring = Array(Int,0)
    for ti in Iterator #TriRing(v,msh)
        face = faces[:,ti]
        cw = (face.==v)[[3,1,2]]
        vi, = face[cw]
        push!(vring,vi)
    end

    p0 = points[:,v]
    vects = Array(Float64,3,length(vring))
    for (j,vi) in enumerate(vring)
        vects[:,j] = points[:,vi] - p0
    end

    #vnormall = vnormal(v,points,faces,Iterator) # If I would like to implement more precise curvature calculation
    R = eye(3)
    C,D,E = ZinchenkoDifertential!(vects,normal,R)

    A = [C D/2;D/2 E]
    k1,k2 = eigvals(-A)
    H = (k1 + k2)/2
    
    return H*2
end

function MeanCurvatures!(curvatures,points,faces,n,VertexVRing)
    for i in 1:size(points,2)
        curvatures[i] = vcurvature(i,points,faces,n[:,i],VertexVRing(i))
    end
end

function FitEllipsoid(points)

    x = points[1,:][:]
    y = points[2,:][:]
    z = points[3,:][:]

    D = [x.*x+y.*y - 2*z.*z x.*x+z.*z-2*y.*y 2*x.*y 2*x.*z 2*y.*z 2*x 2*y 2*z 1+0*x ]

    d2 = x .* x + y .* y + z .* z; #% the RHS of the llsq problem (y's)
    u = ( D' * D ) \ ( D' * d2 );  #% solution to the normal equations

    v = Array(Float64,10)
    v[1] = u[1] +     u[2] - 1;
    v[2] = u[1] - 2 * u[2] - 1;
    v[3] = u[2] - 2 * u[1] - 1;
    v[4:10] = u[3:9];

    # v the 10 parameters describing the ellipsoid / conic algebraically: 
    #               Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
    ### Asuming that ellipsoid is placed at the middle and

    A = v[1]
    B = v[2]
    C = v[3]
    J = v[10]

    a = sqrt(-J/A)
    b = sqrt(-J/B)
    c = sqrt(-J/C)

    return a,b,c
end


function ellipsoid_parameters(points)

    warn("Deprecated syntax use FitEllipsoid instead")
    FitEllipsoid(points)
end

