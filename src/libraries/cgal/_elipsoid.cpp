#include <Python.h>
#include <numpy/arrayobject.h>
#include "elipsoid.h"
#include <iostream>

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

static char module_docstring[] =
    "Here are some scripts which uses CGAL library";
static char test_elipsoid_mesh_docstring[] =
    "Prints the elipsoid mesh to stdout";
static char elipsoid_mesh_docstring[] =
  "Returns elipsoid mesh: verticies, triangles"
  "float angular_bound = 30.0"
  "float radius_bound = 0.1"
  "float distance_bound = 0.1"
  "float radius2_bounding = 2.0";


static PyObject *_test_elipsoid_mesh(PyObject *self, PyObject *args);
static PyObject *_elipsoid_mesh(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"test_elipsoid_mesh", _test_elipsoid_mesh, METH_VARARGS, test_elipsoid_mesh_docstring},
    {"elipsoid_mesh", _elipsoid_mesh, METH_VARARGS, elipsoid_mesh_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC PyInit_elipsoid(void)
{
    
    PyObject *module;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "elipsoid",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    /* Load `numpy` functionality. */
    import_array();

    return module;
}

static PyObject *_test_elipsoid_mesh(PyObject *self, PyObject *args) {

  test_elipsoid_mesh();
  return Py_None;
}

static PyObject *_elipsoid_mesh(PyObject *self, PyObject *args) {

  float a,b,c, angular_bound, radius_bound, distance_bound, radius2_bounding;
  if (!PyArg_ParseTuple(args, "fffffff", &a, &b, &c, &angular_bound, &radius_bound, &distance_bound, &radius2_bounding)) return NULL;

  Mesh ret = elipsoid_mesh(a,b,c,angular_bound,radius_bound,distance_bound,radius2_bounding);
  int *tri = ret.triangles;
  int nt = ret.n_triangles;
  float *vert = ret.verticies;
  int nv = ret.n_verticies;


  npy_intp dim[2] = {nt,3}; // assuming similarity with shape atribute
  PyObject *triangles = PyArray_SimpleNewFromData(2,dim, PyArray_INT,tri);

  npy_intp dimv[2] = {nv,3}; 
  PyObject *verticies = PyArray_SimpleNewFromData(2,dimv, PyArray_FLOAT,vert);

  return Py_BuildValue("OO",verticies,triangles);
}



