#include <stdio.h>
#include <math.h>

#include <vector>
#include <vec.h>
#include <surftrack.h>
#include <subdivisionscheme.h>

//void improvemesh(SurfaceMesh msh,ElTopoParameters elparameters,SurfaceMesh * outmsh)
extern "C" void improvemesh
(
 double * inmsh_verticies,
 int inmsh_Nverticies,
 int * inmsh_triangles,
 int inmsh_Ntriangles,
 double * outmsh_verticies,
 int * outmsh_Nverticies,
 int * outmsh_triangles,
 int * outmsh_Ntriangles,
 double m_proximity_epsilon,
 double m_friction_coefficient,
 double m_min_triangle_area,
 double m_improve_collision_epsilon,
 bool m_use_fraction,
 double m_min_edge_length,
 double m_max_edge_length, 
 double m_max_volume_change,
 double m_min_triangle_angle,
 double m_max_triangle_angle,   
 bool m_use_curvature_when_splitting,
 bool m_use_curvature_when_collapsing,
 double m_min_curvature_multiplier,
 double m_max_curvature_multiplier,
 bool m_allow_vertex_movement,
 double m_edge_flip_min_length_change,
 double m_merge_proximity_epsilon,
 bool m_collision_safety,
 bool m_allow_topology_changes,
 bool m_allow_non_manifold,
 bool m_perform_improvement,
 bool m_verbose
 )
{
  std::vector<Vec3d> vs;
  std::vector<double> masses;

  for ( int i = 0; i < inmsh_Nverticies; ++i )
    {
      vs.push_back( Vec3d( inmsh_verticies[3*i], inmsh_verticies[3*i + 1], inmsh_verticies[3*i + 2] ) );
      masses.push_back( 0.5 );      
    }
    
  std::vector<Vec3st> ts;
  for ( int i = 0; i < inmsh_Ntriangles; ++i )
    {
      ts.push_back( Vec3st( inmsh_triangles[3*i], inmsh_triangles[3*i + 1], inmsh_triangles[3*i + 2] ) );
      printf("%d %d %d\n", inmsh_triangles[3*i], inmsh_triangles[3*i + 1], inmsh_triangles[3*i + 2] );
    }

  
  
  SurfTrackInitializationParameters parameters;

  parameters.m_proximity_epsilon = m_proximity_epsilon;
  parameters.m_friction_coefficient = m_friction_coefficient;
  parameters.m_min_triangle_area = m_min_triangle_area;
  parameters.m_improve_collision_epsilon = m_improve_collision_epsilon;
  parameters.m_use_fraction = m_use_fraction;
  parameters.m_min_edge_length = m_min_edge_length;
  parameters.m_max_edge_length = m_max_edge_length;
  parameters.m_max_volume_change = m_max_volume_change;
  parameters.m_min_triangle_angle = m_min_triangle_angle;
  parameters.m_max_triangle_angle = m_max_triangle_angle;
  parameters.m_use_curvature_when_splitting = m_use_curvature_when_splitting;
  parameters.m_use_curvature_when_collapsing = m_use_curvature_when_collapsing;
  parameters.m_min_curvature_multiplier = m_min_curvature_multiplier;
  parameters.m_max_curvature_multiplier = m_max_curvature_multiplier;
  parameters.m_allow_vertex_movement = m_allow_vertex_movement;
  parameters.m_edge_flip_min_length_change = m_edge_flip_min_length_change;
  parameters.m_merge_proximity_epsilon = m_merge_proximity_epsilon;
  parameters.m_collision_safety = m_collision_safety;
  parameters.m_allow_topology_changes = m_allow_topology_changes;
  parameters.m_allow_non_manifold = m_allow_non_manifold;
  parameters.m_perform_improvement = m_perform_improvement;
  parameters.m_subdivision_scheme = new ButterflyScheme();
  
  SurfTrack surface_tracker( vs, ts, masses, parameters );
  surface_tracker.m_verbose = m_verbose;

  printf("Hello again \n");

  
  surface_tracker.improve_mesh();
  surface_tracker.topology_changes();
  surface_tracker.defrag_mesh();

  // Creating usual array objects
  
  *outmsh_Nverticies = surface_tracker.get_num_vertices();
  //double verticies_out[Nverticies_out*3];
  
  for ( int i = 0; i < *outmsh_Nverticies; ++i )
    {
      outmsh_verticies[3*i+0] = surface_tracker.get_position(i)[0];
      outmsh_verticies[3*i+1] = surface_tracker.get_position(i)[1];
      outmsh_verticies[3*i+2] = surface_tracker.get_position(i)[2];
    }

  *outmsh_Ntriangles = surface_tracker.m_mesh.num_triangles();
//  int triangles_out[Ntriangles_out*3];
  for ( int i = 0; i < *outmsh_Ntriangles; ++i )
    {
      const Vec3st& curr_tri = surface_tracker.m_mesh.get_triangle(i); 
      outmsh_triangles[3*i + 0] = curr_tri[0];
      outmsh_triangles[3*i + 1] = curr_tri[1];
      outmsh_triangles[3*i + 2] = curr_tri[2];
    }

  printf("Nverticies=%d  Ntriangles=%d   volume=%f \n",*outmsh_Nverticies,*outmsh_Ntriangles,surface_tracker.get_volume());
  
//  printmesh(102,Nverticies_out,Ntriangles_out,verticies_out,triangles_out);
//  outmsh->verticies = verticies_out;
  //outmsh->Nverticies = Nverticies_out;
//  outmsh->triangles = triangles_out;
//  outmsh->Ntriangles = Ntriangles_out;

}


extern "C" void hello
(
 int a,
 double b,
 bool s
 )
{
  printf("Hello %f\n",b);
}
