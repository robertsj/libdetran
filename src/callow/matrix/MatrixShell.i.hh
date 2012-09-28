//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MatrixShell.i.hh
 * \author robertsj
 * \date   Sep 20, 2012
 * \brief  MatrixShell.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXSHELL_I_HH_
#define callow_MATRIXSHELL_I_HH_

namespace callow
{

//---------------------------------------------------------------------------//

MatrixShell::MatrixShell(void* context)
  : d_context(context)
{
  /* ... */
}

//---------------------------------------------------------------------------//

MatrixShell::MatrixShell(void* context, const int m, const int n)
  : d_context(context)
{
  set_size(m, n);
}

//---------------------------------------------------------------------------//

MatrixShell::~MatrixShell()
{
  if (d_is_ready)
  {
#ifdef CALLOW_ENABLE_PETSC
    // destroy the petsc matrix.  note, since we constructed it
    // using our own arrays, those still need to be deleted.
    MatDestroy(&d_petsc_matrix);
#endif
  }
}

//---------------------------------------------------------------------------//

void MatrixShell::set_operation()
{
#ifdef CALLOW_ENABLE_PETSC
  MatShellSetOperation(d_petsc_matrix, MATOP_MULT,
                       (void(*)(void))shell_multiply_wrapper);
#endif
  d_is_ready = true;
}

//---------------------------------------------------------------------------//
#ifdef CALLOW_ENABLE_PETSC
inline PetscErrorCode shell_multiply_wrapper(Mat A, Vec x, Vec y)
{
  // get the context and cast
  PetscErrorCode ierr;
  void *context;
  ierr = MatShellGetContext(A, &context);
  Assert(!ierr);
  MatrixShell* foo = (MatrixShell*) context;
  // wrap the petsc vectors
  Vector X(x);
  Vector Y(y);
  // call the actual apply operator.
  foo->multiply(X, Y);
  return ierr;
}
#endif


} // end namespace callow


#endif /* callow_MATRIXSHELL_I_HH_ */
