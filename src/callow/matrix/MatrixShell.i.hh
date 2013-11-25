//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MatrixShell.i.hh
 *  @brief MatrixShell inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef callow_MATRIXSHELL_I_HH_
#define callow_MATRIXSHELL_I_HH_

namespace callow
{

//---------------------------------------------------------------------------//
inline void MatrixShell::set_operation()
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

//----------------------------------------------------------------------------//
//              end of MatrixShell.i.hh
//----------------------------------------------------------------------------//
