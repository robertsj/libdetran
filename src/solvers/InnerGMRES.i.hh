//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InnerGMRES.i.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  InnerGMRES inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef INNERGMRES_I_HH_
#define INNERGMRES_I_HH_

// System
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cstring>

namespace detran
{

template <class D>
void InnerGMRES<D>::solve(int g)
{
  using std::cout;
  using std::endl;

  PetscErrorCode ierr;

  // Set the group for this solve.
  d_g = g;

  if (d_print_out > 0) std::cout << "    Starting SI." << std::endl;

  // Setup boundary conditions.  This sets any conditions fixed for the solve.
  d_boundary->set(g);

  // Set the equations.
  d_sweeper->setup_group(g);

  //-------------------------------------------------------------------------//
  // RIGHT HAND SIDE
  //-------------------------------------------------------------------------//

  // Build the fixed source.  This includes external, fission, and in-scatter.
  d_sweepsource->reset();
  d_sweepsource->build_fixed_with_scatter(g);

  // Sweep.  This is equivalent to
  //   D*inv(L)*Q(Omega).
  State::moments_type B(d_mesh->number_cells(), 0.0);
  d_sweeper->sweep(B);

  // Fill the PETSc Vec with B. Remember to replace it.
  ierr = VecPlaceArray(d_B, &B[0]);
  Insist(!ierr, "Error placing PETSc array");

  //-------------------------------------------------------------------------//
  // SOLVE THE TRANSPORT EQUATION
  //-------------------------------------------------------------------------//

  ierr = KSPSolve(d_solver, d_B, d_X);
  Insist(!ierr, "Error in KSPSolve.");

  //-------------------------------------------------------------------------//
  // POSTPROCESS
  //-------------------------------------------------------------------------//

  // Copy the solution into state.
  double *x_a;
  VecGetArray(d_X, &x_a);
  Insist(!ierr, "Error getting PETSc array.");
  double *phi_a = &d_state->phi(g)[0];
  memcpy(phi_a, x_a, d_state->moments_size()*sizeof(double));
  phi_a = NULL;

  // Get the PETSc residual norm.
  double norm_residal;
  ierr = KSPGetResidualNorm(d_solver, &norm_residal);
  Insist(!ierr, "Error getting residual norm.");

  // Get the number of iterations.
  int iteration;
  ierr = KSPGetIterationNumber(d_solver, &iteration);
  Insist(!ierr, "Error getting iteration number.");

  if (d_print_out > 0)
  {
    printf(" GMRES Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, norm_residal, d_sweeper->number_sweeps());
  }

  if (norm_residal > d_tolerance)
    warning(SOLVER_CONVERGENCE, "    InnerGMRES did not converge.");

  // Replace the storage for B and X.
  ierr = VecResetArray(d_B);
  Insist(!ierr, "Error resetting array.");
  ierr = VecRestoreArray(d_X, &x_a);
  Insist(!ierr, "Error restoring array.");

}

template class InnerGMRES<_1D>;
template class InnerGMRES<_2D>;
template class InnerGMRES<_3D>;

} // namespace detran

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

inline PetscErrorCode apply_WGTO_1D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as InnerGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::InnerGMRES<detran::_1D> *inner =
    (detran::InnerGMRES<detran::_1D>*) ctx;
  // Call the actual apply operator.
  return inner->apply_WGTO(A, x, y);
}

inline PetscErrorCode apply_WGTO_2D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as InnerGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::InnerGMRES<detran::_2D> *inner =
    (detran::InnerGMRES<detran::_2D>*) ctx;
  // Call the actual apply operator.
  return inner->apply_WGTO(A, x, y);
}

inline PetscErrorCode apply_WGTO_3D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as InnerGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::InnerGMRES<detran::_3D> *inner =
    (detran::InnerGMRES<detran::_3D>*) ctx;
  // Call the actual apply operator.
  return inner->apply_WGTO(A, x, y);
}

#endif /* INNERGMRES_I_HH_ */

//---------------------------------------------------------------------------//
//              end of InnerGMRES.i.hh
//---------------------------------------------------------------------------//
