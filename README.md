# SurfaceGeometry

Library for triangular surface mesh generation, stabilisation and surface parameter estimation. It is intended to be used as surface velocity integrator for boundary element method, so I it also includes vertex normal, curvature calculation and also circulators. 

## Future directions

- [ ] Implenent general kind of mesh generator with CGAL by passing signed distance function from julia 
- [ ] Sphere mesh generation with Octahedron
- [ ] Mesh subdivision method
- [ ] Interfacing remesher from http://www.gris.informatik.tu-darmstadt.de/~sfuhrman/remesher.html which would solve bad initial mesh form CGAL
- [ ] Interface for Eltopo library allowing to handle more extreme deformations

## Mesh Generation

For now only ellipsoid mesh can be generated with Distmesh or CGAL library. Distmesh accepts one parameter `step` which charecterizes the edge length of triangular element which is stored in object with type `DistmeshSurfaceMesher`. The interface to genrate the mesh so is
```
a,b,c = 2,1,1
mesher = DistmeshSurfaceMesher()
# And if you want to use CGAL insted
# mesher = CGALSurfaceMesher()
points, faces = EllipsoidMesh(a,b,c,mesher)
```
Similarly you can use CGAL library only changing 2nd line of code above. CGAL mesher has four parameters - `AngularBound` the minimum allowed angle (Warning: stability for CGAL mesher is being (proven?) observed below 30 degrees), `RadiusBound`, `DistanceBound` measures the distance from faces to the meshed surface, `BoundingRadius` sets the bounds where object is going to be meshed. The interface to set theese parameters is `CGALSurfaceMesher(AngularBound,RadiusBound,DistanceBound,BoundingRadius)` or `CGALSurfaceMesher(;AngularBound=25)`. 

## Estimation of surface parameters

For calculation of vertex normals and mean curvature functions `NormalVectors!` and `MeanCurvatures` are created. A basic interface to them is
```
n = Array(size(points)...)
NormalVectors!(n,points,faces,i->FaceVRing(i,faces))

curvatures = Array(size(points,2))
MeanCurvatures!(curvatures,points,faces,n,i->FaceVRing(i,faces))
```
where in the arguments `VertexVRing` and `FaceVRring` are vertex and face circulators discoussed more in next section. 

For surface global parameters use function `volume(points,faces)`, for prolate ellipsoid fit of the points use `FitEllipsoid(points)`.

## Data structure and circulators

Circulators are iterators which read mesh topology (for example `faces`) and ansver the questions such as which verticies are directly connected to given vertex, which traingular faces has common vertex. The vertex circulator `VertexVRing`, face circulator `FaceVRing` and oposite edge circulator `DoubleVertexVRing` can be used as follows
```
vertex = 1
for i in VertexVRing(vertex,faces)
    println("This is vertex $i")
end

for i in FaceVRing(vertex,faces)
    println("This is face $i which has common vertex $vertex")
end

for (v1,v2) in DoubleVertexVring(vertex,faces)
    println("The oposite edge to $vertex in counterclockwise direction is $v1 $v2")
end 
```

The circulators above uses seach over face elements, which in some cases are too slow. Instead you can save all searches beforehand in data structure. For example you can initialise face based data structure and use all previuos circulator functions
```
dataStructure = FaceBasedDS(faces)

vertex = 1
for i in VertexVRing(vertex,dataStructure)
    println("This is vertex $i")
end
```
Or you can choose connectivity data structure which for me have prooven to be by a facor of 10 faster
```
valence = 7
dataStructure = ConnectivityDS(faces,valence)
```

## Surface velocity field integration and stabilisation

Assuming that velocity function is given as `velocity(t,points,faces)` the integration for points with Euler method would look like
```
t = 0
h = 0.1
p = points

for i in 1:10
    v = velocity(t,p,faces)
    p += h*v
end 
```
In this way it is however more likely to experience the degeneration of triangles as surface evolves. To deal with it we can seperate normal and tangential velocity where the later one can be adjusted for keeping good mesh.

Here we use "kinetic" energy minimization function to find tangential velocity as [Zinchenko1997] and [Zinchenko2013] had used (the springs with damping also are quite popular [Cristini]). The interface to integrate above with stabilisation is as follows
```
t = 0
h = 0.1
p = points

for i in 1:10
    v = velocity(t,p,faces)

    res =  stabilise(points,faces,v)
    println("Energy before minimization Finit=$(res.Finit) after Fres=$(res.Fres)")
    v = res.vres
    
    p += h*v
end 
```

Eugene (incompressable velocity field proposed in [Eigene]]) test is being used for checking a need for stabilisation. For example without stabilisation we have a following picture:

| ![](https://raw.githubusercontent.com/akels/SurfaceGeometry.jl/master/img/OriginalField.gif) | ![](https://raw.githubusercontent.com/akels/SurfaceGeometry.jl/master/img/NormalField.gif) | ![](https://github.com/akels/SurfaceGeometry.jl/blob/master/img/StabilisedField.gif) |
|---|---|---|
| | | |

