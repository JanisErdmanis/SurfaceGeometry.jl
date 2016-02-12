#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This is the file where you can look for an example of iterating, geting basic vertex positions, outputing triangles
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef typename Tr::Vertex_handle Vertex_handle;
typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
typedef typename Tr::Point Point;
typedef typename Tr::Facet Facet;
typedef typename Tr::Edge Edge;

float a,b,c;

FT elipsoid (Point_3 p) {

  const FT x = p.x();
  const FT y = p.y();
  const FT z = p.z();
  
  //const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x*x/a/a + y*y/b/b + z*z/c/c -1;
}

struct Mesh
{
  float verticies[100000];
  int triangles[100000];
  int n_verticies;
  int n_triangles;
};

Mesh mesh_repr(const C2t3& c2t3) {
  const Tr& tr = c2t3.triangulation();

  Mesh r;
  std::map<Vertex_handle, int> V;
  int inum = 0;
  Finite_vertices_iterator vit = tr.finite_vertices_begin();
  while(vit != tr.finite_vertices_end()) {

    // making an integer representation of vertex pointers
    V[vit] = inum;

    // obtaining vertex positions from vertex pointer vit
    Point p = static_cast<Point>(vit->point());
    r.verticies[inum*3 + 0] = p.x();
    r.verticies[inum*3 + 1] = p.y();
    r.verticies[inum*3 + 2] = p.z();
      
    ++vit;
    inum+=1;
  }
  r.n_verticies = inum; 

  // Orienting surface
  bool success = true;
  Finite_facets_iterator fit = tr.finite_facets_begin();
  std::set<Facet> oriented_set;
  std::stack<Facet> stack;

  typename Tr::size_type number_of_facets = c2t3.number_of_facets();
  
  //  CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

  while (oriented_set.size() != number_of_facets) 
    {
      while ( fit->first->is_facet_on_surface(fit->second) == false ||
              oriented_set.find(*fit) != oriented_set.end() ||
		
              oriented_set.find(c2t3.opposite_facet(*fit)) !=
              oriented_set.end() ) 
        {
          ++fit;
        }
      oriented_set.insert(*fit);
      stack.push(*fit);
      while(! stack.empty() )
        {
          Facet f = stack.top();
          stack.pop();
          for(int ih = 0 ; ih < 3 ; ++ih) {
            const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
            const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));

            const typename C2t3::Face_status face_status
              = c2t3.face_status(Edge(f.first, i1, i2));
            if(face_status == C2t3::REGULAR) {
              Facet fn = c2t3.neighbor(f, ih);
              if (oriented_set.find(fn) == oriented_set.end()) {
                if(oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
                  {
                    oriented_set.insert(fn);
                    stack.push(fn);
                  }
                else {
                  success = false; // non-orientable
                }
              }
            }
            else if(face_status != C2t3::BOUNDARY) {
              success = false; // non manifold, thus non-orientable
            }
          } // end "for each neighbor of f"
        } // end "stack non empty"
    } // end "oriented_set not full"
   
  // Outputting it [orineted_set]

  
  //Finite_facets_iterator fit = tr.finite_facets_begin();
  
  inum = 0;
  for (typename std::set<Facet>::const_iterator fit = oriented_set.begin();fit != oriented_set.end(); ++fit) {

    typename Tr::Cell_handle cell = fit->first;
    const int& index = fit->second;

    r.triangles[inum*3 + 0] = V[cell->vertex(tr.vertex_triple_index(index, 0))];
    r.triangles[inum*3 + 1] = V[cell->vertex(tr.vertex_triple_index(index, 1))];
    r.triangles[inum*3 + 2] = V[cell->vertex(tr.vertex_triple_index(index, 2))];
    
//    ++fit;
    inum +=1;
  }
  r.n_triangles = inum;
  
  return r;
}

Mesh elipsoid_mesh(float aa, float bb, float cc, float angular_bound, float radius_bound, float distance_bound, float radius2_bounding) {

  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  a = aa;
  b = bb;
  c = cc;
  
  Surface_3 surface(elipsoid, Sphere_3(CGAL::ORIGIN, radius2_bounding)); 
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angular_bound,radius_bound,distance_bound);
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  Mesh result = mesh_repr(c2t3);
  return result;
}

int test_elipsoid_mesh() {

  float angular_bound = 30.0;
  float radius_bound = 0.1;
  float distance_bound = 0.1;
  float radius2_bounding = 2.0;
  float a = 2;
  float b = 1;
  float c = 1;
  
  Mesh result = elipsoid_mesh(a,b,c,angular_bound,radius_bound,distance_bound,radius2_bounding);

  int i = 0;
  while (i != result.n_triangles) {
    std::cout << result.triangles[i*3 + 0] << " " << result.triangles[i*3 + 1] << " " << result.triangles[i*3 + 2] << std::endl;
    i+=1;
  }

  i = 0;
  while (i != result.n_verticies) {
    std::cout << result.verticies[i*3 + 0] << " " << result.verticies[i*3 + 1] << " " << result.verticies[i*3 + 2] << std::endl;
    i+=1;
  }
  
}





