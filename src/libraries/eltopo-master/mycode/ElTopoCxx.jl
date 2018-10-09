#g++ ../eltopo3d/obj/*.o -o eltopo.so -fPIC -shared -llapack -lblas -lstdc++ -lm -I../common -I../eltopo3d

module ElTopo

using Parameters

@with_kw struct SurfTrackInitializationParameters
    # Elements closer than this are considered "near" (or proximate)
    proximity_epsilon::Float64 = 1e-4
    friction_coefficient::Float64 = 0.0
    min_triangle_area::Float64 = 1e-7

    # Collision epsilon to use during mesh improvment operations (i.e. if any mesh elements are closer than this, the operation is 
    # aborted).  NOTE: This should be greater than collision_epsilon, to prevent improvement operations from moving elements into 
    # a collision configuration.
    improve_collision_epsilon::Float64 = 2e-6
    
    # Whether to set the min and max edge lengths as fractions of the initial average edge length
    use_fraction::Bool = false
    
    # If use_fraction is true, these are taken to be fractions of the average edge length of the new surface.
    # If use_fraction is false, these are absolute.
    min_edge_length::Float64 = 0.05
    max_edge_length::Float64 = 0.2 
    max_volume_change::Float64 = 0.1
    
    # In-triangle angles to enforce
    min_triangle_angle::Float64 = 0
    max_triangle_angle::Float64 = 180   
    
    use_curvature_when_splitting::Bool = false
    use_curvature_when_collapsing::Bool = false
    
    # Clamp curvature scaling to these values
    min_curvature_multiplier::Float64 = 1
    max_curvature_multiplier::Float64 = 1
    
    allow_vertex_movement::Bool = false
    
    # Minimum edge length improvement in order to flip an edge
    edge_flip_min_length_change::Float64 = 0.05
    
    # Elements within this distance will trigger a merge attempt   
    merge_proximity_epsilon::Float64 = 1e-5
    
    # Whether to enforce collision-free surfaces (including during mesh maintenance operations)
    collision_safety::Bool = true
    
    # Whether to allow changes in topology
    allow_topology_changes::Bool = false
    
    # Wether to allow non-manifold (edges incident on more than two triangles)
    allow_non_manifold::Bool = false
    
    # Whether to allow mesh improvement
    perform_improvement::Bool = true
end


using Cxx
Libdl.dlopen(pwd()*"/eltopo.so",Libdl.RTLD_GLOBAL)
addHeaderDir("../eltopo3d",kind=C_System)
addHeaderDir("../common",kind=C_System)

cxxinclude("vector")
cxxinclude("subdivisionscheme.h")
cxxinclude("surftrack.h")

############# C++ code ############

### Seems that moving Array to C++ and backwards is not so easally achievable due to 
### Its usage of stack allocated arrays for points. Something similar as StaicArrays.
### Therefore I need manually to manually upload vertices to the mesh for execution. 

### After that it is possible to
### icxx"$parameters.m_proximity_epsilon;"
function constructparameters(p::SurfTrackInitializationParameters)

    icxx"""
      SurfTrackInitializationParameters parameters;
      parameters.m_proximity_epsilon = $(p.proximity_epsilon);
      parameters.m_friction_coefficient = $(p.friction_coefficient);
      parameters.m_min_triangle_area = $(p.min_triangle_area);
      parameters.m_improve_collision_epsilon = $(p.improve_collision_epsilon);
      parameters.m_use_fraction = $(p.use_fraction);
      parameters.m_min_edge_length = $(p.min_edge_length);
      parameters.m_max_edge_length = $(p.max_edge_length);
      parameters.m_max_volume_change = $(p.max_volume_change);
      parameters.m_min_triangle_angle = $(p.min_triangle_angle);
      parameters.m_max_triangle_angle = $(p.max_triangle_angle);
      parameters.m_use_curvature_when_splitting = $(p.use_curvature_when_splitting);
      parameters.m_use_curvature_when_collapsing = $(p.use_curvature_when_collapsing);
      parameters.m_min_curvature_multiplier = $(p.min_curvature_multiplier);
      parameters.m_max_curvature_multiplier = $(p.max_curvature_multiplier);
      parameters.m_allow_vertex_movement = $(p.allow_vertex_movement);
      parameters.m_edge_flip_min_length_change = $(p.edge_flip_min_length_change);
      parameters.m_merge_proximity_epsilon = $(p.merge_proximity_epsilon);
      parameters.m_collision_safety = $(p.collision_safety);
      parameters.m_allow_topology_changes = $(p.allow_topology_changes);
      parameters.m_allow_non_manifold = $(p.allow_non_manifold);
      parameters.m_perform_improvement = $(p.perform_improvement);
      parameters.m_subdivision_scheme = new ButterflyScheme();
    parameters;
    """
    # icxx"""
    #   SurfTrackInitializationParameters parameters;
    #   parameters.m_proximity_epsilon = 1e-4;
    #   parameters.m_friction_coefficient = 0.0;
    #   parameters.m_min_triangle_area = 1e-7;
    #   parameters.m_improve_collision_epsilon = 2e-6;
    #   parameters.m_use_fraction = false;
    #   parameters.m_min_edge_length = 0.05;
    #   parameters.m_max_edge_length = 0.2;
    #   parameters.m_max_volume_change = 0.1;
    #   parameters.m_min_triangle_angle = 35;
    #   parameters.m_max_triangle_angle = 180;
    #   parameters.m_use_curvature_when_splitting = false;
    #   parameters.m_use_curvature_when_collapsing = false;
    #   parameters.m_min_curvature_multiplier = 1.0;
    #   parameters.m_max_curvature_multiplier = 1.0;
    #   parameters.m_allow_vertex_movement = false;
    #   parameters.m_edge_flip_min_length_change = 1e-8;
    #   parameters.m_merge_proximity_epsilon = 1e-5;
    #   parameters.m_collision_safety = true;
    #   parameters.m_allow_topology_changes = false;
    #   parameters.m_allow_non_manifold = false;
    #   parameters.m_perform_improvement = true;
    #   parameters.m_subdivision_scheme = new ButterflyScheme();
    # parameters;
    # """
end

function improve_mesh(points,faces,parameters::SurfTrackInitializationParameters)

    ### Here one actually would do conversion of parameters
    p = constructparameters(parameters)
    faces = faces .- 1
    
    ### Creating C++ vector objects to construct C++ object
    
    vs = icxx"std::vector<Vec<3,double>> vs; vs;"
    ms = icxx"std::vector<double> masses; masses;"
    ts = icxx"std::vector<Vec<3,size_t>> ts; ts;"

    for i in 1:size(points,2)
        p1,p2,p3 = points[1,i], points[2,i], points[3,i]
        icxx"$vs.push_back(Vec3d($p1,$p2,$p3));"
        icxx"$ms.push_back(0.5);"
    end

    for i in 1:size(faces,2)
        v1,v2,v3 = faces[1,i], faces[2,i], faces[3,i]
        icxx"$ts.push_back(Vec3st($v1,$v2,$v3));"
    end


    ### Defining outputs
    vsout = icxx"std::vector<Vec<3,double>> vsout; vsout;"
    tsout = icxx"std::vector<Vec<3,size_t>> tsout; tsout;"

    ### Constructing surface tracker and creating outputs
    icxx"""
    SurfTrack surface_tracker($vs,$ts,$ms,$p);
    surface_tracker.m_verbose = false;

    surface_tracker.improve_mesh();
    surface_tracker.topology_changes();
    surface_tracker.defrag_mesh();

    for ( int i = 0; i < surface_tracker.get_num_vertices(); ++i ) 
            $vsout.push_back(surface_tracker.get_position(i));

    for ( int i = 0; i < surface_tracker.m_mesh.num_triangles(); ++i ) 
            $tsout.push_back(surface_tracker.m_mesh.get_triangle(i));

    """

    ### Construction of a propper julia vectors
    ### I will use map this time!!!

    vsoutj = map(x->Float64(icxx"$vsout[$(x[1])][$(x[2])];"),Iterators.product(0:(Int(icxx"$vsout.size();")-1),0:2))
    tsoutj = map(x->Int(UInt64(icxx"$tsout[$(x[1])][$(x[2])];")),Iterators.product(0:(Int(icxx"$tsout.size();")-1),0:2))
    
    vsoutj,tsoutj .+ 1

end

export SurfTrackInitializationParameters, improve_mesh

end


using Main.ElTopo
using JLD
@load "sphere.jld"

p = SurfTrackInitializationParameters()





