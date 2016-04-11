module SurfaceGeometry

include("MeshGenerationMethods/Utils.jl")
include("MeshGenerationMethods/CGAL.jl")
include("MeshGenerationMethods/Distmesh.jl")
export EllipsoidMesh, SurfaceMesh, volume, DistmeshSurfaceMesher, CGALSurfaceMesher, subdivision

include("Iterators.jl")
include("ComplexDS.jl")
export FaceBasedDS, ConnectivityDS
export FaceVRing, VertexVRing, DoubleVertexVRing

include("Properties.jl")
export NormalVectors!, MeanCurvatures!, FitEllipsoid, volume

### Essential from stabilistation
include("StabilisationMethods/stabilisationV1.jl")
include("StabilisationMethods/stabilisationV2.jl")
include("StabilisationMethods/stabilisationV3.jl")
include("StabilisationMethods/Stabilisation.jl")
include("StabilisationMethods/ElTopo.jl")

export stabilise, improvemesh, Elparameters, improvemeshcol

end # module
