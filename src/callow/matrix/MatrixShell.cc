//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MatrixShell.cc
 *  @brief MatrixShell member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts

 */
//----------------------------------------------------------------------------//

#include "MatrixShell.hh"

namespace callow
{

//----------------------------------------------------------------------------//
MatrixShell::MatrixShell(void* context)
  : d_context(context)
{
  /* ... */
}

//----------------------------------------------------------------------------//
MatrixShell::MatrixShell(void* context, const int m, const int n)
  : d_context(context)
{
  set_size(m, n);
}

//----------------------------------------------------------------------------//
MatrixShell::~MatrixShell()
{
  if (d_is_ready)
  {
//#ifdef CALLOW_ENABLE_PETSC
//    // destroy the petsc matrix.  note, since we constructed it
//    // using our own arrays, those still need to be deleted.
//    MatDestroy(&d_petsc_matrix);
//#endif
  }
}


} // end namespace callow




