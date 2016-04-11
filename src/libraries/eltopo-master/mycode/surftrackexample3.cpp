// g++ surftrackexample3.cpp libeltopo_release.a -llapack -lblas -lstdc++ -lm -I../common -I../eltopo3d -o  test 
//C hello world example
#include <stdio.h>
//#include <stdbool.h>
#include <math.h>

#include <vector>
#include <vec.h>
#include <surftrack.h>
#include <subdivisionscheme.h>

#define NSIZE 20000

#include "mesh.c"

int printmesh(int ii, int Nverticies, int Ntriangles, double * verticies, int * triangles) {

  char buf[256];
  snprintf(buf, sizeof buf, "%s%d%s", "output/", ii, ".obj");

  FILE *fp;

  fp = fopen(buf, "w+");

  //fprintf(fp, "This is testing for fprintf...\n");
  

  for (int i=0;i<Nverticies;i++){
    fprintf(fp,"v %f %f %f \n",verticies[3*i],verticies[3*i+1],verticies[3*i+2]);
  }

  for (int i=0;i<Ntriangles;i++) {
    fprintf(fp,"f %d %d %d\n",triangles[3*i]+1,triangles[3*i+1]+1,triangles[3*i+2]+1);
  }
    

  //  fputs("This is testing for fputs...\n", fp);
  fclose(fp);
}

struct ElTopoParameters
{
    ElTopoParameters()
    {
      m_collision_safety = true;
      m_allow_non_manifold = false;
      m_allow_topology_changes = false;
      m_allow_vertex_movement = false;
      m_use_fraction = false;
      m_verbose = false;
      m_use_curvature_when_splitting = false;
      m_use_curvature_when_collapsing = false;
      m_perform_improvement = true;

      m_proximity_epsilon = 1e-4;
      m_friction_coefficient = 0.0;
      m_min_triangle_area = 1e-7;
      m_improve_collision_epsilon = 2e-6;
      m_min_curvature_multiplier = 1;
      m_max_curvature_multiplier = 1;
      m_edge_flip_min_length_change = 1e-8;
      m_merge_proximity_epsilon = 1e-5;
        
      m_max_volume_change = 0.1;
      m_min_triangle_angle = 35;
      m_max_triangle_angle = 180;
      m_min_edge_length = 0.01;
      m_max_edge_length = 0.005;
    }

  /// Elements closer than this are considered "near" (or proximate)
  double m_proximity_epsilon;
    
  double m_friction_coefficient;
    
  double m_min_triangle_area;
    
  /// Collision epsilon to use during mesh improvment operations (i.e. if any mesh elements are closer than this, the operation is 
  /// aborted).  NOTE: This should be greater than collision_epsilon, to prevent improvement operations from moving elements into 
  /// a collision configuration.
  double m_improve_collision_epsilon;
    
  /// Whether to set the min and max edge lengths as fractions of the initial average edge length
  bool m_use_fraction;
    
  /// If use_fraction is true, these are taken to be fractions of the average edge length of the new surface.
  /// If use_fraction is false, these are absolute.
  double m_min_edge_length;
  double m_max_edge_length; 
  double m_max_volume_change;
    
  // In-triangle angles to enforce
  double m_min_triangle_angle;
  double m_max_triangle_angle;   
    
    
  bool m_use_curvature_when_splitting;
  bool m_use_curvature_when_collapsing;
    
  // Clamp curvature scaling to these values
  double m_min_curvature_multiplier;
  double m_max_curvature_multiplier;
    
  bool m_allow_vertex_movement;
    
  /// Minimum edge length improvement in order to flip an edge
  double m_edge_flip_min_length_change;
    
  /// Elements within this distance will trigger a merge attempt   
  double m_merge_proximity_epsilon;
    
  /// Whether to enforce collision-free surfaces (including during mesh maintenance operations)
  bool m_collision_safety;
    
  /// Whether to allow changes in topology
  bool m_allow_topology_changes;
    
  /// Wether to allow non-manifold (edges incident on more than two triangles)
  bool m_allow_non_manifold;
    
  /// Whether to allow mesh improvement
  bool m_perform_improvement;

  // Printing out a bunch of stuff
  bool m_verbose;
};

struct SurfaceMesh
{
  double * verticies;
  int Nverticies;
  int * triangles;
  int Ntriangles;
};

void improvemesh(SurfaceMesh msh,ElTopoParameters elparameters,SurfaceMesh * outmsh)
{
  std::vector<Vec3d> vs;
  std::vector<double> masses;
    
  for ( int i = 0; i < msh.Nverticies; ++i )
    {
      vs.push_back( Vec3d( msh.verticies[3*i], msh.verticies[3*i + 1], msh.verticies[3*i + 2] ) );
      masses.push_back( 0.5 );      
    }
    
  std::vector<Vec3st> ts;
  for ( int i = 0; i < msh.Ntriangles; ++i )
    {
      ts.push_back( Vec3st( msh.triangles[3*i], msh.triangles[3*i + 1], msh.triangles[3*i + 2] ) );
    }
  
  SurfTrackInitializationParameters parameters;

  parameters.m_proximity_epsilon = elparameters.m_proximity_epsilon;
  parameters.m_friction_coefficient = elparameters.m_friction_coefficient;
  parameters.m_min_triangle_area = elparameters.m_min_triangle_area;
  parameters.m_improve_collision_epsilon = elparameters.m_improve_collision_epsilon;
  parameters.m_use_fraction = elparameters.m_use_fraction;
  parameters.m_min_edge_length = elparameters.m_min_edge_length;
  parameters.m_max_edge_length = elparameters.m_max_edge_length;
  parameters.m_max_volume_change = elparameters.m_max_volume_change;
  parameters.m_min_triangle_angle = elparameters.m_min_triangle_angle;
  parameters.m_max_triangle_angle = elparameters.m_max_triangle_angle;
  parameters.m_use_curvature_when_splitting = elparameters.m_use_curvature_when_splitting;
  parameters.m_use_curvature_when_collapsing = elparameters.m_use_curvature_when_collapsing;
  parameters.m_min_curvature_multiplier = elparameters.m_min_curvature_multiplier;
  parameters.m_max_curvature_multiplier = elparameters.m_max_curvature_multiplier;
  parameters.m_allow_vertex_movement = elparameters.m_allow_vertex_movement;
  parameters.m_edge_flip_min_length_change = elparameters.m_edge_flip_min_length_change;
  parameters.m_merge_proximity_epsilon = elparameters.m_merge_proximity_epsilon;
  parameters.m_collision_safety = elparameters.m_collision_safety;
  parameters.m_allow_topology_changes = elparameters.m_allow_topology_changes;
  parameters.m_allow_non_manifold = elparameters.m_allow_non_manifold;
  parameters.m_perform_improvement = elparameters.m_perform_improvement;

  parameters.m_subdivision_scheme = new ButterflyScheme();
  
  SurfTrack surface_tracker( vs, ts, masses, parameters );
  surface_tracker.m_verbose = elparameters.m_verbose;

  surface_tracker.improve_mesh();
  surface_tracker.topology_changes();
  surface_tracker.defrag_mesh();

  // Creating usual array objects
  
  int Nverticies_out = surface_tracker.get_num_vertices();
  //double verticies_out[Nverticies_out*3];
  
  for ( int i = 0; i < Nverticies_out; ++i )
    {
      outmsh->verticies[3*i+0] = surface_tracker.get_position(i)[0];
      outmsh->verticies[3*i+1] = surface_tracker.get_position(i)[1];
      outmsh->verticies[3*i+2] = surface_tracker.get_position(i)[2];
    }

  int Ntriangles_out = surface_tracker.m_mesh.num_triangles();
//  int triangles_out[Ntriangles_out*3];
  for ( int i = 0; i < Ntriangles_out; ++i )
    {
      const Vec3st& curr_tri = surface_tracker.m_mesh.get_triangle(i); 
      outmsh->triangles[3*i + 0] = curr_tri[0];
      outmsh->triangles[3*i + 1] = curr_tri[1];
      outmsh->triangles[3*i + 2] = curr_tri[2];
    }

  printf("Nverticies=%d  Ntriangles=%d   volume=%f \n",Nverticies_out,Ntriangles_out,surface_tracker.get_volume());
  
//  printmesh(102,Nverticies_out,Ntriangles_out,verticies_out,triangles_out);

//  outmsh->verticies = verticies_out;
  outmsh->Nverticies = Nverticies_out;
//  outmsh->triangles = triangles_out;
  outmsh->Ntriangles = Ntriangles_out;

}

int main()
{

  // Learn to write interface with cppwrapper
  // Interface this code
  // Put this in seperate library ElTopo.jl
    
  printmesh(0,Nverticies,Ntriangles,verticies,triangles);
   
  SurfaceMesh msh;
  msh.verticies = verticies;
  msh.Nverticies = Nverticies;
  msh.triangles = triangles;
  msh.Ntriangles = Ntriangles;

  ElTopoParameters elparameters = ElTopoParameters();

  double verticiesout[NSIZE];
  int trianglesout[NSIZE];
  SurfaceMesh mshout;
  mshout.verticies = verticiesout;
  mshout.triangles = trianglesout;
  
  improvemesh(msh,elparameters,&mshout);

  printmesh(101,mshout.Nverticies,mshout.Ntriangles,mshout.verticies,mshout.triangles);
 
  //printf("%f\n",elparameters.m_max_volume_change);
  
  printf("Hello world\n");
  return 0;
}
