//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Operator.cc
 * \brief  Operator member definitions
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#include "Operator.hh"
#include "utilities/DBC.hh"

namespace detran
{

Operator::Operator(const size_t m,
                   const size_t n)
  : d_number_rows(m)
  , d_number_columns(n)
  , d_is_assembled(false)
{
  // Preconditions
  Require(m > 0);
  Require(n > 0);

  // Set the matrix size.  This applies to all operators, which should
  // know their sizes from the start.
  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_SELF, &d_A);
  ierr = MatSetSizes(d_A, m, n, PETSC_DETERMINE, PETSC_DETERMINE);

  // Postconditions
  Require(!ierr);
}

Operator::~Operator()
{
  MatDestroy(&d_A);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Operator.cc
//---------------------------------------------------------------------------//
