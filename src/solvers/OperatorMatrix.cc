//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   OperatorMatrix.cc
 * \brief  OperatorMatrix 
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#include "OperatorMatrix.hh"

namespace detran
{

OperatorMatrix::OperatorMatrix(const size_t m, const size_t n)
{

  // Create appropriate matrix type.
  PetscErrorCode ierr;
  ierr = MatSetType(d_A, MATSEQAIJ);

  // Build
  build();

  // Postconditions
  Ensure(!ierr);
}

void OperatorMatrix::assemble()
{
  if (!d_is_assembled)
  {
    MatAssemblyBegin(d_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(d_A, MAT_FINAL_ASSEMBLY);
    d_is_assembled = true;
  }
}

// Default implementation.  Note, shell matrices have to do something else!
void OperatorMatrix::display(const int output, const std::string name) const
{
  if (output == STDOUT)
  {
    MatView(d_A, PETSC_VIEWER_STDOUT_SELF);
    return;
  }
  PetscViewer viewer;
  if (output == BINARY)
    PetscViewerBinaryOpen(PETSC_COMM_SELF, name.c_str(), FILE_MODE_WRITE, &viewer);
  else if (output == ASCII)
    PetscViewerASCIIOpen(PETSC_COMM_SELF, name.c_str(), &viewer);
  else
    THROW("Invalid output switch for Matrix display");
  MatView(d_A, viewer);
  PetscViewerDestroy(&viewer);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file OperatorMatrix.cc
//---------------------------------------------------------------------------//
