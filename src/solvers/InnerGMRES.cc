//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerGMRES.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerGMRES class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// Detran
#include "InnerGMRES.hh"
//
//#include "CMR.hh"

namespace detran
{

// Constructor
template <class D>
InnerGMRES<D>::InnerGMRES(SP_input          input,
                          SP_state          state,
                          SP_mesh           mesh,
                          SP_material       material,
                          SP_quadrature     quadrature,
                          SP_boundary       boundary,
                          SP_externalsource q_e,
                          SP_fissionsource  q_f)
  :  InnerIteration<D>::InnerIteration(input,
                                       state,
                                       mesh,
                                       material,
                                       quadrature,
                                       boundary,
                                       q_e,
                                       q_f)
{

  // Error from PETSc calls.
  PetscErrorCode ierr;

  // Determine the sizes of the moments.
  int d_moments_size = mesh->number_cells();

  // Determine the sizes of any reflected boundary fluxes.
  int d_boundary_size = 0*quadrature->number_angles()/2;

  int number_unknowns = d_moments_size + d_boundary_size;

  ierr = MatCreateShell(PETSC_COMM_WORLD,
                        number_unkowns,
                        number_unkowns,
                        number_unkowns,
                        number_unkowns,
                        this,
                        &d_operator);
  Insist(!ierr, "Error creating MR shell matrix.");

  // Create the corresponding vectors.
  ierr = VecCreateSeq(PETSC_COMM_WORLD,
                      number_unkowns,
                      &d_X);
  ierr = VecCreateSeq(PETSC_COMM_WORLD,
                      number_unkowns,
                      &d_B);

  // Set the operator.
  ierr = set_operation();
  Insist(!ierr, "Error setting matrix-vector operation.");

  // Create the KSP object.
  ierr = KSPCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating KSP object.");

  // Set the operator.
  KSPSetOperators(d_solver, d_operator, d_operator, SAME_NONZERO_PATTERN);

  // Set tolerances.
  ierr = KSPSetTolerances(d_solver,
                          d_tolerance,
                          PETSC_DEFAULT,
                          PETSC_DEFAULT,
                          PETSC_DEFAULT);

  // Allow for command line flags.
  ierr = KSPSetFromOptions(d_solver);

}

template <class D>
PetscErrorCode InnerGMRES<D>::apply_WGTO(Mat A, Vec x, Vec y)
{

  PetscErrorCode ierr;

  // Get the array from the Krylov vector x.
  double *x_a;
  ierr = VecGetArray(x, &x_a); CHKERRQ(ierr);

  // Assign the array to a std::vector.  Is there a memory-efficient
  // way to do this?
  int n = d_state->moments_size();
  State::moments_type x_v(n);
  for (int i = 0; i < n; i++)
  {
    x_v[i] = x_a[i];
  }

  // y <-- x
  State::moments_type y_v(x_v);

  // Build the within group source.  This is equivalent to
  //   x <-- M*S*x
  d_sweepsource->reset();
  d_sweepsource->build_within_group_scatter(d_g, x_v);

  // Set the incident boundary fluxes.

  // Sweep over space and angle.  This is equivalent to
  //   x <-- D*inv(L)*M*S*x
  d_sweeper->sweep(x_v);

  // Return the following:
  //  y <-- x - D*inv(L)*M*S*x = (I-D*inv(L)*M*S)*x <-- A*x
  for (int i = 0; i < n; i++)
  {
    y_v[i] -= x_v[i];
  }
  double *y_a;
  ierr = VecGetArray(y, &y_a); CHKERRQ(ierr);
  for (int i = 0; i < n; i++)
  {
    y_a[i] = y_v[i];
  }

  // Restore the arrays.
  ierr = VecRestoreArray(x, &x_a); CHKERRQ(ierr);
  ierr = VecRestoreArray(y, &y_a); CHKERRQ(ierr);

  // No error.
  return 0;
}

template <class D>
PetscErrorCode InnerGMRES<D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_WGTO_3D);
}
template <>
PetscErrorCode InnerGMRES<_2D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_WGTO_2D);
}
template <>
PetscErrorCode InnerGMRES<_1D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_WGTO_1D);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of InnerGMRES.cc
//---------------------------------------------------------------------------//
