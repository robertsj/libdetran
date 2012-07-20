//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSLEPc.cc
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  EigenSLEPc class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#include "detran_config.h"

#ifdef DETRAN_ENABLE_SLEPC

// Detran
#include "EigenSLEPc.hh"

// System
#include <iostream>

namespace detran
{

// Constructor
template <class D>
EigenSLEPc<D>::EigenSLEPc(SP_input          input,
                          SP_state          state,
                          SP_mesh           mesh,
                          SP_material       material,
                          SP_quadrature     quadrature,
                          SP_boundary       boundary,
                          SP_fissionsource  q_f)
  : Base(input, state, mesh, material, quadrature, boundary, q_f)
  , d_size(mesh->number_cells())
  , d_mg_solves(0)
{

  PetscErrorCode ierr;

  // Create the shell matrix.
  ierr = MatCreateShell(PETSC_COMM_WORLD, d_size, d_size, d_size, d_size,
                        this, &d_operator);
  Insist(!ierr, "Error creating shell matrix.");

  // Set the operator.
  ierr = set_operation();
  Insist(!ierr, "Error setting matrix-vector operation.");

  // Create the context.
  ierr = EPSCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating EPS context.");

  // Set the operator.
  ierr = EPSSetOperators(d_solver, d_operator, PETSC_NULL);
  Insist(!ierr, "Error setting EPS operator.");

  // Set the problem type.
  ierr = EPSSetProblemType(d_solver, EPS_NHEP);
  Insist(!ierr, "Error setting EPS problem type.");

  // Set the solver type to krylovschur and one largest eigenvalue.
  ierr = EPSSetType(d_solver, EPSKRYLOVSCHUR);
  Insist(!ierr, "Error defaulting EPS to EPSKRYLOVSCHUR.");
  ierr = EPSSetWhichEigenpairs(d_solver, EPS_LARGEST_MAGNITUDE);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_solver);
  Insist(!ierr, "Error setting EPS from options.");

}

template <class D>
PetscErrorCode EigenSLEPc<D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_eigen_3D);
}
template <>
PetscErrorCode EigenSLEPc<_2D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_eigen_2D);
}
template <>
PetscErrorCode EigenSLEPc<_1D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_eigen_1D);
}

} // end namespace detran

#endif

//---------------------------------------------------------------------------//
//              end of EigenSLEPc.cc
//---------------------------------------------------------------------------//
