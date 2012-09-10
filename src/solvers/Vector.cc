//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.cc
 * \brief  Vector member definitions
 * \author Jeremy Roberts
 * \date   Sep 8, 2012
 */
//---------------------------------------------------------------------------//

#include "Vector.hh"

namespace detran
{

Vector::Vector(const size_t m, const double val)
  : d_size(m)
  , d_is_assembled(false)
{
  // Preconditions
  Require(m > 0);

  // Create the vector
  PetscErrorCode ierr;
  ierr = VecCreateSeq(PETSC_COMM_SELF, m, &d_V);
  ierr = VecSet(d_V, val);

  // Get acces to and then "restore" the array, but not actually.
  // PETSc requires the restore call to ensure
  // safe coding, but by passing null, we get to keep it.
  // We'll code safely...
  ierr = VecGetArray(d_V, &d_array);
  ierr = VecRestoreArray(d_V, PETSC_NULL);

  // Because we initialize to zero, we can assemble right
  // away, since we might just be filling this from a multiplication.
  assemble();

  // Postconditions
  Ensure(!ierr);
}

Vector::~Vector()
{
  VecDestroy(&d_V);
}

void Vector::assemble()
{
  if (!d_is_assembled)
  {
    VecAssemblyBegin(d_V);
    VecAssemblyEnd(d_V);
    d_is_assembled = true;
  }
}

void Vector::display() const
{
  VecView(d_V, PETSC_VIEWER_STDOUT_SELF);
}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Vector.cc
//---------------------------------------------------------------------------//
