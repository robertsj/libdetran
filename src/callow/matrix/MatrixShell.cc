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
  /* ... */
}

//----------------------------------------------------------------------------//
void MatrixShell::set_size(const int m, const int n)
{
  Require(m > 0);
  d_m = m;
  d_n = n;
  if (!d_n) d_n = d_m;
  d_sizes_set = true;
#ifdef DETRAN_ENABLE_PETSC
  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_SELF, &d_petsc_matrix);
  ierr = MatSetSizes(d_petsc_matrix, d_m, d_n, PETSC_DETERMINE, PETSC_DETERMINE);
  ierr = MatSetType(d_petsc_matrix, MATSHELL);
  ierr = MatShellSetContext(d_petsc_matrix, d_context);
  ierr = MatSetUp(d_petsc_matrix);
  Ensure(!ierr);
#endif
  set_operation();
  d_is_ready = true;
}

//----------------------------------------------------------------------------//
void MatrixShell::display(bool forceprint) const
{
  std::cout << "MatrixShell:" << std::endl
            << "  # rows = " << d_m << std::endl
            << "  # cols = " << d_n << std::endl
            << std::endl;
}

} // end namespace callow


//----------------------------------------------------------------------------//
//              end of MatrixShell.cc
//----------------------------------------------------------------------------//
