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
inline void InnerGMRES<D>::solve(int g)
{
  using std::cout;
  using std::endl;


  PetscErrorCode ierr;

  // Set the group for this solve.
  d_g = g;

  if (d_print_out > 0) std::cout << "    Starting GMRES." << std::endl;

  // Set the preconditioner, if used.
  if (d_use_pc) d_pc->set_group(g);

  // Set the equations.
  d_sweeper->setup_group(g);

  //-------------------------------------------------------------------------//
  // BUILD RIGHT HAND SIDE
  //-------------------------------------------------------------------------//

  // Right hand side vector.
  State::moments_type B(d_mesh->number_cells(), 0.0);

  // Build the right hand side.
  build_rhs(B);

  //-------------------------------------------------------------------------//
  // SOLVE THE TRANSPORT EQUATION
  //-------------------------------------------------------------------------//

  ierr = KSPSolve(d_solver, d_B, d_X);
  Insist(!ierr, "Error in KSPSolve.");

  //-------------------------------------------------------------------------//
  // POSTPROCESS
  //-------------------------------------------------------------------------//

  // Copy the flux solution into state.
  double *X_a;
  VecGetArray(d_X, &X_a);
  Insist(!ierr, "Error getting PETSc array.");
  double *phi_a = &d_state->phi(g)[0];
  memcpy(phi_a, X_a, d_moments_size*sizeof(double));

//  if (d_boundary->has_reflective())
//  {
//    // Set the incident boundary flux.
//    d_boundary->set_incident(g, X_a + d_moments_size);
//  }

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
  //ierr = VecResetArray(d_B);
  //Insist(!ierr, "Error resetting array.");
  ierr = VecRestoreArray(d_X, &X_a);
  Insist(!ierr, "Error restoring array.");

}

template <class D>
inline void InnerGMRES<D>::build_rhs(State::moments_type &B)
{

  // Clear the group boundary.  This is because the right hand side should
  // be based only on fixed boundaries.
  d_boundary->clear(d_g);

  // Setup boundary conditions.  This sets any conditions fixed for the solve.
  d_boundary->set(d_g);

  // Build the fixed source.  This includes external, fission, and in-scatter.
  d_sweepsource->reset();
  d_sweepsource->build_fixed_with_scatter(d_g);

  // If no reflective, one sweep and out.
  if (!d_boundary->has_reflective())
  {
    d_sweeper->sweep(B);
  }
  // Otherwise, solve the reflecting condition problem.
  else
  {
    // Create dummy moment container for convergence.
    State::moments_type B_old(B.size(), 0.0);
    // Set sweeper to update boundaries on the fly.
    d_sweeper->set_update_boundary(true);
    int iteration;
    double error = 0;
    for (iteration = 0; iteration < 1090; iteration++)
    {
      // Update boundary.  This updates boundaries due to reflection, etc.
      d_boundary->update(d_g);
      // Swap old and new.
      std::swap(B, B_old);
      // Sweep.
      d_sweeper->sweep(B);
      // Flux residual using L-infinity.
      error = norm_residual(B_old, B, "Linf");
      if (error < d_tolerance)
      {
        d_reflective_solve_iterations = iteration + 1;
        break;
      }
    }
    d_sweeper->set_update_boundary(false);
  }

  // Fill the PETSc Vec with B. Remember to replace it.
  double *B_a;
  VecSet(d_B, 0.0);
  VecGetArray(d_B, &B_a);
  for (int i = 0; i < d_moments_size; i++)
  {
    B_a[i] = B[i];
  }
  VecRestoreArray(d_B, &B_a);

}

//---------------------------------------------------------------------------//
// SPECIALIZED OPERATORS
//---------------------------------------------------------------------------//

template <class D>
inline PetscErrorCode InnerGMRES<D>::apply_WGTO(Mat A, Vec X, Vec Y)
{

  using std::cout;
  using std::endl;
  PetscErrorCode ierr;

  // Get the array from the Krylov vector X.
  double *X_a;
  ierr = VecGetArray(X, &X_a); CHKERRQ(ierr);

  // Assing array to phi_vac.
  State::moments_type phi_original(d_moments_size, 0.0);
  State::moments_type phi_update(d_moments_size, 0.0);

  for (int i = 0; i < d_moments_size; i++)
  {
    phi_original[i] = X_a[i];
    phi_update[i]   = X_a[i];
  }

  // Reset the source and place the original outgoing boundary flux.
  d_boundary->clear(d_g);

  if (d_boundary->has_reflective())
  {
    // *** SET THE INCIDENT BOUNDARY ***
    d_boundary->set_incident(d_g, X_a + d_moments_size);
  }

  // Reset the source to zero.
  d_sweepsource->reset();
  d_sweepsource->build_within_group_scatter(d_g, phi_original);
  // Sweep.  This gives X <-- D*inv(L)*M*S*X
  d_sweeper->sweep(phi_update);

  // ******** BUILD THE OUTGOING VECTOR

  // Get the array for the outgoing vector.
  double *Y_a;
  ierr = VecGetArray(Y, &Y_a); CHKERRQ(ierr);

  // Assign the moment values.  This gives X <- (I-D*inv(L)*M*S)*X
  for (int i = 0; i < d_moments_size; i++)
  {
    Y_a[i] = phi_original[i] - phi_update[i];
  }

  if (d_boundary->has_reflective())
  {
    // Update the boundary and fetch.
    d_boundary->update(d_g);
    vec_dbl psi_update(d_boundary_size, 0.0);

    // *** EXTRACT THE INCIDENT BOUNDARY ***
    d_boundary->get_incident(d_g, &psi_update[0]);

    // Add the boundary values.
    for (int a = 0; a < d_boundary_size; a++)
    {
      ///  face, o, a, d_g  = psi_out;
      Y_a[a + d_moments_size] =
        X_a[a + d_moments_size] - psi_update[a];
    }
  }
  // Restore the arrays.
  ierr = VecRestoreArray(X, &X_a); CHKERRQ(ierr);
  ierr = VecRestoreArray(Y, &Y_a); CHKERRQ(ierr);

  // No error.
  return 0;
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
