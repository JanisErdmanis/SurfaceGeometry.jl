module SurfaceGeometry

using LinearAlgebra
eye(n) = Matrix(1.0I, n, n)

include("Iterators.jl")
include("ComplexDS.jl")
export FaceBasedDS, ConnectivityDS
export FaceVRing, VertexVRing, DoubleVertexVRing

include("Utils.jl")
export subdivision, volume

include("Properties.jl")
export NormalVectors!, MeanCurvatures!, FitEllipsoid

### Essential from stabilistation
include("StabilisationMethods/stabilisationV1.jl")
include("StabilisationMethods/stabilisationV2.jl")
include("StabilisationMethods/stabilisationV3.jl")
include("StabilisationMethods/Stabilisation.jl")
export stabilise

### Just redirecting
using ElTopo
export Elparameters, improvemesh, improvemeshcol

end # module










