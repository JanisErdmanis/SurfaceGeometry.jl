#include("Utils.jl")

immutable CGALSurfaceMesher
    AngularBound::Float64
    RadiusBound::Float64
    DistanceBound::Float64
    BoundingRadius::Float64
end

### For an easy syntax
function CGALSurfaceMesher(;AngularBound=20,RadiusBound=0.5,DistanceBound=0.1,BoundingRadius=8.0)
    CGALSurfaceMesher(AngularBound,RadiusBound,DistanceBound,BoundingRadius)
end

function EllipsoidMesh(a,b,c,cg::CGALSurfaceMesher)

    fdis(x,y,z) = x^2/a^2 + y^2/b^2 + z^2/c^2 - 1
    SurfaceMesh(fdis,cg)
end

function SurfaceMesh(fdis::Function,ds::CGALSurfaceMesher)
    
    const sdf_c = cfunction(fdis, Cfloat, (Cfloat,Cfloat,Cfloat))

    verticies = zeros(Float32,10000)
    faces = zeros(Int32,10000)

    Nverticies = Ref{Cint}(0)
    Nfaces = Ref{Cint}(0)

    cgallib = Pkg.dir("SurfaceGeometry","src","libraries","cgal","libcgalmesh.so")

    ### This part is really ugly
    eval(:(ccall((:genmesh, $cgallib),Void,(Ptr{Void},Cfloat,Cfloat,Cfloat,Cfloat,Ptr{Cfloat},Ptr{Cint},Ref{Cint},Ref{Cint}),$sdf_c,$(ds.AngularBound),$(ds.RadiusBound),$(ds.DistanceBound),$(ds.BoundingRadius),$verticies,$faces,$Nverticies,$Nfaces)))
    
    #ccall((:genmesh, "/home/janiserdmanis/Documents/cgal/libelipsoid.so"),Void,(Ptr{Void},Cfloat,Cfloat,Cfloat,Cfloat,Ptr{Cfloat},Ptr{Cint},Ref{Cint},Ref{Cint}),sdf_c,ds.AngularBound,ds.RadiusBound,ds.DistanceBound,ds.BoundingRadius,verticies,faces,Nverticies,Nfaces)

    verticies = reshape(verticies[1:(3*Nverticies[])],3,Int(Nverticies[]))
    faces = reshape(faces[1:(3*Nfaces[])],3,Int(Nfaces[])) + 1


    ### This part generally can be solved with volume integration
    if !isoriented(verticies,faces)
        faces[1,:], faces[2,:] = faces[2,:], faces[1,:]
    end

    verticies,faces = map(Float64,verticies), map(Int,faces)

    for vi in 1:size(verticies,2)
        verticies[:,vi] = pushback(x -> fdis(x[1],x[2],x[3]),verticies[:,vi])
    end

    return verticies,faces
end

