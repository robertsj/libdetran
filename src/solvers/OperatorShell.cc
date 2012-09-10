//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OperatorShell.cc
 * \brief  OperatorShell 
 * \author Jeremy Roberts
 * \date   Sep 9, 2012
 */
//---------------------------------------------------------------------------//

#include "OperatorShell.hh"

namespace detran
{

OperatorShell::OperatorShell(const size_t m,
                             const size_t n,
                             void* context)
  : Operator(m, n)
{
  // No preconditions

  // Set matrix type and context
  PetscErrorCode ierr;
  ierr = MatSetType(d_A, MATSHELL);
  ierr = MatShellSetContext(d_A, context);

  // There is no assembly for a shell.
  d_is_assembled = true;

  // Set the matrix shell operations.
  ierr = set_operation();

  // Postconditions
  Ensure(!ierr);
}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file OperatorShell.cc
//---------------------------------------------------------------------------//
