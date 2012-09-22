// allowing callow to be its own python module
%module pycallow
%{

#include "callow_config.hh"
//
#include "utils/SP.hh"
#include "callow/utils/DBC.hh"
#include "callow/utils/Initialization.hh"
#include "callow/utils/Typedefs.hh"
//
#include "callow/vector/Vector.hh"
//
#include "callow/matrix/MatrixBase.hh"
#include "callow/matrix/Matrix.hh"
#include "callow/matrix/MatrixShell.hh"
//
#include "callow/solver/LinearSolverCreator.hh"
#include "callow/solver/Richardson.hh"
#include "callow/solver/Jacobi.hh"
#include "callow/solver/GaussSeidel.hh"
#include "callow/solver/GMRES.hh"
#ifdef CALLOW_ENABLE_PETSC
#include "callow/solver/PetscSolver.hh"
#endif
//
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/preconditioner/PCJacobi.hh"
#include "callow/preconditioner/PCILU0.hh"
#include "callow/preconditioner/PCShell.hh"

%} // end module pycallow

%include "callow.i"

// Anyhere in C/C++ that we need (argc, argv), 
%typemap(in) (int argc, char *argv[]) 
{
  /* Check if is a list */
  if (PyList_Check($input)) 
  {
    int i;
    $1 = PyList_Size($input);
    $2 = (char **) malloc(($1+1)*sizeof(char *));
    for (i = 0; i < $1; i++) 
    {
      PyObject *o = PyList_GetItem($input,i);
      if (PyString_Check(o))
        $2[i] = PyString_AsString(PyList_GetItem($input,i));
      else 
      {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($2);
        return NULL;
      }
    }
    $2[i] = 0;
  } 
  else 
  {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap(freearg) (int argc, char *argv[]) {
  free((char *) $2);
}
