// g++ bareinterface.cpp libeltopo_release.a -llapack -lblas -lm -I../common -I../eltopo3d -o  test 
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#define NSIZE 20000

#include <eltopo.h>
#include <subdivisionscheme.h>

struct ElTopoMesh input; 
struct ElTopoMesh output; 
struct ElTopoGeneralOptions general_otions;
struct ElTopoStaticOperationsOptions options; 
struct ElTopoDefragInformation defrag_info;

#include "mesh.c"

int main()
{

  general_otions.m_verbose = 1;
  general_otions.m_collision_safety = 1;
  general_otions.m_proximity_epsilon = 1e-3; //!!!!!!!!!

  options.m_perform_improvement = true;
  options.m_allow_topology_changes = false;
  options.m_max_volume_change = 0.1;
  options.m_min_edge_length = 0.5;
  options.m_max_edge_length = 1.5;
  options.m_min_triangle_area = 1;
  options.m_min_triangle_angle = 20;
  options.m_max_triangle_angle = 170;
  options.m_use_curvature_when_splitting = true;
  options.m_use_curvature_when_collapsing = true;
  options.m_min_curvature_multiplier = 1.0;
  options.m_max_curvature_multiplier = 1.0;
  options.m_allow_vertex_movement = true;
  options.m_edge_flip_min_length_change = 0.001;
  options.m_merge_proximity_epsilon = 0.0001;
  options.m_subdivision_scheme = new ButterflyScheme();
  options.m_collision_safety = true;
  options.m_allow_non_manifold = false;

  int vertex_is_remove[NSIZE];
  int vertex_index[NSIZE];
  int split_edge[2*NSIZE];
  int triangle_is_remove[NSIZE];
  int triangle_index[NSIZE];
  int new_tri[3*NSIZE];
  int defragged_triangle_map[NSIZE];
  int defragged_vertex_map[NSIZE];
  defrag_info.num_vertex_changes = 0;
  defrag_info.vertex_is_remove = vertex_is_remove;
  defrag_info.vertex_index = vertex_index;
  defrag_info.split_edge = split_edge;
  defrag_info.num_triangle_changes = 0;
  defrag_info.triangle_is_remove = triangle_is_remove;
  defrag_info.triangle_index = triangle_index;
  defrag_info.new_tri = new_tri;
  defrag_info.defragged_triangle_map_size = 0;
  defrag_info.defragged_triangle_map = defragged_triangle_map;
  defrag_info.defragged_vertex_map_size = 0;
  defrag_info.defragged_vertex_map = defragged_vertex_map;

  double vertex_masses[Nverticies];
  for(int i=0;i<Nverticies;i++) {
    vertex_masses[i] = 1;
  }
  
  input.num_vertices = Nverticies;
  input.vertex_locations = verticies;
  input.num_triangles = Ntriangles;
  input.triangles = triangles;
  input.vertex_masses = vertex_masses;

  double out_vertices[NSIZE];
  int out_triangles[NSIZE];
  double out_masses[NSIZE];
  output.vertex_locations = out_vertices;
  output.triangles = out_triangles;
  output.vertex_masses = out_masses;

  // Testing out interface

  el_topo_static_operations(&input,&general_otions,&options,&defrag_info,&output);

  printf("Hello world\n");
  return 0;
}
