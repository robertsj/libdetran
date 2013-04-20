//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   callow.i
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

%import "detran_utilities.i"
%include "callow_config.hh"

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

//---------------------------------------------------------------------------//
// definitions
//---------------------------------------------------------------------------//

%include "utils/CallowDefinitions.hh"

//---------------------------------------------------------------------------//
// vector
//---------------------------------------------------------------------------//

%include "vector/Vector.i"

//---------------------------------------------------------------------------//
// matrix
//---------------------------------------------------------------------------//

%include "matrix/Matrix.i"

//---------------------------------------------------------------------------//
// linear solver
//---------------------------------------------------------------------------//

%include "preconditioner/Preconditioner.i"
%include "solver/Solver.i"

