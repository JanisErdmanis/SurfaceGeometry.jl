maxtheta = 0

for ti in 1:size(f,2)
    v1,v2,v3 = f[:,ti]
    pc = (p[:,v1] + p[:,v2] + p[:,v3])/3
    nx = cross(p[:,v2]-p[:,v1],p[:,v3]-p[:,v1])

    
    theta = acos(dot(nx,pc)/norm(nx)/norm(pc))
    if theta>maxtheta
        maxtheta = theta
    end

    println(theta)
end

### Checking if calculated normals are directed to the outside

info("Checking a normal vector calculation")

f = map(Int,f)


normals = Array(Float64,size(p)...);
NormalVectors!(normals,p,f,i->FaceVRing(i,f))    

maxtheta = 0

for ti in 1:size(f,2)
    v1,v2,v3 = f[:,ti]
    pc = (p[:,v1] + p[:,v2] + p[:,v3])/3
    nx = (normals[:,v1] + normals[:,v2] + normals[:,v3])/3
    
    theta = acos(dot(nx,pc)/norm(nx)/norm(pc))
    if theta>maxtheta
        maxtheta = theta
    end

    println(theta)
end
