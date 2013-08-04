//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_callow.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for callow library.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") callow
%{
#include <stddef.h>
#include "callow/callow_config.hh"
#include "callow/utils/Initialization.hh"
#include "callow/utils/Typedefs.hh"
#include "callow/vector/Vector.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/matrix/MatrixShell.hh"
#include "callow/matrix/MatrixDense.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/EigenSolverCreator.hh"
#include "callow/preconditioner/Preconditioner.hh"
%}

// Hide templates from SWIG
%inline
{
#define CALLOW_EXPORT
#define CALLOW_TEMPLATE_EXPORT(...)
#define CALLOW_INSTANTIATE_EXPORT(...)
}

%import "detran_utilities.i"
//%include "callow_config.hh"

//---------------------------------------------------------------------------//
// initialization
//---------------------------------------------------------------------------//

void callow_initialize(int argc, char *argv[]);
void callow_finalize();

//---------------------------------------------------------------------------//
// setup for numerical arrays
//---------------------------------------------------------------------------//

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
  import_array();
%}
%numpy_typemaps(double, NPY_DOUBLE, int)
%numpy_typemaps(in,     NPY_INT, int)
%apply (int* IN_ARRAY1, double* IN_ARRAY2, int DIM1) 
       {(int *j, double *v, int n)}
%apply (int* IN_ARRAY1, int* IN_ARRAY2, double* IN_ARRAY3, int DIM1) 
       {(int *i, int *j, double* v, int n)}
// Vector views.  These give us direct access to C arrays of 
// basic types via Numpy arrays.  We'll apply these below to 
// expose C++ std::vectors of those types.
%apply (double** ARGOUTVIEW_ARRAY1, int *DIM1) 
       {(double** a, int *n)}
%apply (int** ARGOUTVIEW_ARRAY1, int *DIM1) 
       {(int** a, int *n)}

//---------------------------------------------------------------------------//
// definitions
//---------------------------------------------------------------------------//

//%include "utils/CallowDefinitions.hh"

//---------------------------------------------------------------------------//
// vector
//---------------------------------------------------------------------------//

//%include "vector/Vector.i"

//---------------------------------------------------------------------------//
// matrix
//---------------------------------------------------------------------------//

//%include "matrix/Matrix.i"

//---------------------------------------------------------------------------//
// linear solver
//---------------------------------------------------------------------------//

//%include "preconditioner/Preconditioner.i"
//%include "solver/Solver.i"


%inline
{
  // Vectors to arrays (and then to Numpy). Careful, or you'll end up
  // in memory management pergatory.  These are used for plotting, etc.
  // Usage:
  //   v = vec_dbl(10, 1.23)
  //   a = vec_asarray(v)
  //   a[0] = 2.34
  //   # do something in Numpy with a.  Then for good measure, do
  //   del a
  //   # and no matter what, don't use a after v is out of scope!
  void vec_asarray(std::vector<double> &v, double **a, int *n)
  {
    *n = v.size();
    *a = &v[0];
  }
  void vec_asarray(std::vector<int> &v, int **a, int *n)
  {
    *n = v.size();
    *a = &v[0];
  }
}