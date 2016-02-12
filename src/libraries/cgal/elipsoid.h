struct Mesh
{
  float verticies[100000];
  int triangles[100000];
  int n_verticies;
  int n_triangles;
};

int test_elipsoid_mesh();
Mesh elipsoid_mesh(float a, float b, float c, float angular_bound, float radius_bound, float distance_bound, float radius2_bounding);
