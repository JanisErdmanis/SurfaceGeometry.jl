module SurfaceGeometry

include("MeshGeneration.jl")
export pushback, ellipsoid_mesh_matlab, ellipsoid_mesh_cgal

include("Iterators.jl")
export FaceRing, find_triangle_vertex, find_other_triangle_edge

include("Properties.jl")
export vnormal, ZinchenkoDifertential!, vcurvature, ellipsoid_parameters, volume

include("Integrator.jl")
export TimeStepers, EilerStep, AdamsStep, step!, properstep

include("ViewMesh.jl")
export view3d
# package code goes here

end # module
