NSIZE = 100000

using Parameters
@with_kw type Elparameters
    #/// Elements closer than this are considered "near" (or proximate)
    m_proximity_epsilon::Float64 = 1e-4
    m_friction_coefficient::Float64 = 0.0
    m_min_triangle_area::Float64 = 1e-7

    #/// Collision epsilon to use during mesh improvment operations (i.e. if any mesh elements are closer than this, the operation is 
    #/// aborted).  NOTE: This should be greater than collision_epsilon, to prevent improvement operations from moving elements into 
    #/// a collision configuration.
    m_improve_collision_epsilon::Float64 = 2e-6
    
    #/// Whether to set the min and max edge lengths as fractions of the initial average edge length
    m_use_fraction::Bool = false
    
    #/// If use_fraction is true, these are taken to be fractions of the average edge length of the new surface.
    #/// If use_fraction is false, these are absolute.
    m_min_edge_length::Float64 = 0.05
    m_max_edge_length::Float64 = 0.2 
    m_max_volume_change::Float64 = 0.1
    
    #// In-triangle angles to enforce
    m_min_triangle_angle::Float64 = 0
    m_max_triangle_angle::Float64 = 180   
    
    m_use_curvature_when_splitting::Bool = false
    m_use_curvature_when_collapsing::Bool = false
    
    #// Clamp curvature scaling to these values
    m_min_curvature_multiplier::Float64 = 1
    m_max_curvature_multiplier::Float64 = 1
    
    m_allow_vertex_movement::Bool = false
    
    #/// Minimum edge length improvement in order to flip an edge
    m_edge_flip_min_length_change::Float64 = 0.05
    
    #/// Elements within this distance will trigger a merge attempt   
    m_merge_proximity_epsilon::Float64 = 1e-5
    
    #/// Whether to enforce collision-free surfaces (including during mesh maintenance operations)
    m_collision_safety::Bool = true
    
    #/// Whether to allow changes in topology
    m_allow_topology_changes::Bool = false
    
    #/// Wether to allow non-manifold (edges incident on more than two triangles)
    m_allow_non_manifold::Bool = false
    
    #/// Whether to allow mesh improvement
    m_perform_improvement::Bool = true

    #// Printing out a bunch of stuff
    m_verbose::Bool = false

    m_dt::Float64 = 0
end

const libpath = Pkg.dir("SurfaceGeometry","src","libraries","eltopo-wrapper") * "/eltopo.so"

function improvemesh(verticies,triangles,par)
    triangles = map(Int32,triangles) - 1
    
    outmsh_verticies = zeros(Float64,NSIZE)
    outmsh_triangles = zeros(Int32,NSIZE)
    outmsh_Nverticies = Ref{Cint}(0)
    outmsh_Ntriangles = Ref{Cint}(0)

    # @show ccall((:improvemesh, "/home/janiserdmanis/Documents/eltopo-master-old/mycode/eltopo.so"),Void,
    #             (Ptr{Cdouble},Cint,Ptr{Cint},Cint,
    #              Ptr{Cdouble},Ref{Cint},Ptr{Cint},Ref{Cint}),
    #             inmsh_verticies,length(inmsh_verticies),inmsh_triangles,length(inmsh_triangles),
    # outmsh_verticies,outmsh_Nverticies,outmsh_triangles,outmsh_Ntriangles)

    
    ccall((:improvemesh, libpath),Void,
                (Ptr{Cdouble},Cint,Ptr{Cint},Cint,
                 Ptr{Cdouble},Ref{Cint},Ptr{Cint},Ref{Cint},
                 Cdouble,Cdouble,Cdouble,Cdouble,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cint,Cint,Cdouble,Cdouble,Cint,Cdouble,Cdouble,Cint,Cint,Cint,Cint,Cint
                 ),
                verticies,size(verticies,2),triangles,size(triangles,2),
    outmsh_verticies,outmsh_Nverticies,outmsh_triangles,outmsh_Ntriangles,
    par.m_proximity_epsilon,
    par.m_friction_coefficient,
    par.m_min_triangle_area,
    par.m_improve_collision_epsilon,
    par.m_use_fraction,
    par.m_min_edge_length,
    par.m_max_edge_length,
    par.m_max_volume_change,
    par.m_min_triangle_angle,
    par.m_max_triangle_angle,
    par.m_use_curvature_when_splitting,
    par.m_use_curvature_when_collapsing,
    par.m_min_curvature_multiplier,
    par.m_max_curvature_multiplier,
    par.m_allow_vertex_movement,
    par.m_edge_flip_min_length_change,
    par.m_merge_proximity_epsilon,
    par.m_collision_safety,
    par.m_allow_topology_changes,
    par.m_allow_non_manifold,
    par.m_perform_improvement,
    par.m_verbose
    )
    
    pout = outmsh_verticies[1:3*outmsh_Nverticies[]]
    pout = reshape(pout,3,convert(Int,outmsh_Nverticies[]))
    tout = outmsh_triangles[1:3*outmsh_Ntriangles[]]
    tout = reshape(tout,3,convert(Int,outmsh_Ntriangles[]))

    return pout,map(Int,tout)+1
end


function improvemeshcol(verticies,triangles,newverticies,par)
    triangles = map(Int32,triangles) - 1
    actual_dt = Ref{Float64}(0)
    
    outmsh_verticies = zeros(Float64,NSIZE)
    outmsh_triangles = zeros(Int32,NSIZE)
    outmsh_Nverticies = Ref{Cint}(0)
    outmsh_Ntriangles = Ref{Cint}(0)

    # @show ccall((:improvemesh, "/home/janiserdmanis/Documents/eltopo-master-old/mycode/eltopo.so"),Void,
    #             (Ptr{Cdouble},Cint,Ptr{Cint},Cint,
    #              Ptr{Cdouble},Ref{Cint},Ptr{Cint},Ref{Cint}),
    #             inmsh_verticies,length(inmsh_verticies),inmsh_triangles,length(inmsh_triangles),
    # outmsh_verticies,outmsh_Nverticies,outmsh_triangles,outmsh_Ntriangles)

    
    ccall((:improvecol, libpath),Void,
                (Ptr{Cdouble},Cint,Ptr{Cint},Cint,
                 Ptr{Cdouble},Ref{Cint},Ptr{Cint},Ref{Cint},
                 Cdouble, Ptr{Cdouble}, Ref{Cdouble},
                 Cdouble,Cdouble,Cdouble,Cdouble,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cint,Cint,Cdouble,Cdouble,Cint,Cdouble,Cdouble,Cint,Cint,Cint,Cint,Cint
                 ),
                verticies,size(verticies,2),triangles,size(triangles,2),
    outmsh_verticies,outmsh_Nverticies,outmsh_triangles,outmsh_Ntriangles,
    par.m_dt, newverticies, actual_dt,
    par.m_proximity_epsilon,
    par.m_friction_coefficient,
    par.m_min_triangle_area,
    par.m_improve_collision_epsilon,
    par.m_use_fraction,
    par.m_min_edge_length,
    par.m_max_edge_length,
    par.m_max_volume_change,
    par.m_min_triangle_angle,
    par.m_max_triangle_angle,
    par.m_use_curvature_when_splitting,
    par.m_use_curvature_when_collapsing,
    par.m_min_curvature_multiplier,
    par.m_max_curvature_multiplier,
    par.m_allow_vertex_movement,
    par.m_edge_flip_min_length_change,
    par.m_merge_proximity_epsilon,
    par.m_collision_safety,
    par.m_allow_topology_changes,
    par.m_allow_non_manifold,
    par.m_perform_improvement,
    par.m_verbose
    )
    
    pout = outmsh_verticies[1:3*outmsh_Nverticies[]]
    pout = reshape(pout,3,convert(Int,outmsh_Nverticies[]))
    tout = outmsh_triangles[1:3*outmsh_Ntriangles[]]
    tout = reshape(tout,3,convert(Int,outmsh_Ntriangles[]))

    return actual_dt[],pout,map(Int,tout)+1
end

