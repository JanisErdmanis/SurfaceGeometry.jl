module SurfaceGeometry

include("MeshGenerationMethods/Utils.jl")
include("MeshGenerationMethods/CGAL.jl")
include("MeshGenerationMethods/Distmesh.jl")
export EllipsoidMesh, SurfaceMesh, volume, DistmeshSurfaceMesher, CGALSurfaceMesher

include("Iterators.jl")
include("ComplexDS.jl")
export FaceBasedDS, ConnectivityDS
export FaceVRing, VertexVRing, DoubleVertexVRing

include("Properties.jl")
export NormalVectors!, MeanCurvatures!, FitEllipsoid,

i### Essential from stabilistation
include("StabilisationMethods/stabilisationV1.jl")
include("StabilisationMethods/stabilisationV2.jl")
include("StabilisationMethods/Stabilisation.jl")
export stabilise

end # module
