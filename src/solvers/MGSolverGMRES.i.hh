//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MGSolverGMRES.i.hh
 * \author robertsj
 * \date   Jun 19, 2012
 * \brief  MGSolverGMRES inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGSOLVERGMRES_I_HH_
#define detran_MGSOLVERGMRES_I_HH_

#include "utilities/MathUtilities.hh"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cstring>

namespace detran
{

// Public Interface

template <class D>
inline void MGSolverGMRES<D>::solve(const double keff)
{

  using std::cout;
  using std::endl;

  PetscErrorCode ierr;

//  cout << " MGSolverGMRES solver... " << endl;
//  cout << " upscatter_cutoff = " << d_upscatter_cutoff << endl;
//  cout << "    number_groups = " << d_number_groups << endl;
//  cout << "   upscatter_size = " << d_upscatter_size << endl;
//  cout << "     moments_size = " << d_moments_size << endl;
//  cout << "    boundary_size = " << d_boundary_size << endl;
//  cout << "   moments_size_g = " << d_moments_size_group << endl;
//  cout << "  boundary_size_g = " << d_boundary_size_group << endl;

  double norm_resid = 0.0;
  int iteration = 0;

  //-------------------------------------------------------------------------//
  // DOWNSCATTER BLOCK
  //-------------------------------------------------------------------------//

  for (int g = 0; g < d_upscatter_cutoff; g++)
  {
    d_wg_solver->solve(g);
  }
  // Is there anything to say about the residual?

  //-------------------------------------------------------------------------//
  // UPSCATTER BLOCK
  //-------------------------------------------------------------------------//

  if (d_upscatter_size > 0)
  {

    //-----------------------------------------------------------------------//
    // BUILD RIGHT HAND SIDE
    //-----------------------------------------------------------------------//
    State::vec_moments_type B(d_upscatter_size,
                              State::moments_type(d_moments_size_group, 0.0));
    build_rhs(B);

    //-----------------------------------------------------------------------//
    // SOLVE MULTIGROUP TRANSPORT EQUATION
    //-----------------------------------------------------------------------//

    ierr = KSPSolve(d_solver, d_B, d_X);
    Insist(!ierr, "Error in KSPSolve.");

    //-----------------------------------------------------------------------//
    // POST-PROCESS
    //-----------------------------------------------------------------------//

    // Copy the flux solution into state.
    double *X_a;
    VecGetArray(d_X, &X_a);
    Insist(!ierr, "Error getting PETSc array.");
    for (int g = d_upscatter_cutoff; g < d_number_groups; g++)
    {
      double *phi_a = &d_state->phi(g)[0];
      int offset = (g - d_upscatter_cutoff) *
                   (d_moments_size_group + d_boundary_size_group);
      memcpy(phi_a, X_a + offset, d_moments_size_group*sizeof(double));
      phi_a = NULL;
    }

    // Replace the storage for X.
    ierr = VecRestoreArray(d_X, &X_a);
    Insist(!ierr, "Error restoring array.");

    // Get the PETSc residual norm.
    ierr = KSPGetResidualNorm(d_solver, &norm_resid);
    Insist(!ierr, "Error getting residual norm.");

    // Get the number of iterations.
    ierr = KSPGetIterationNumber(d_solver, &iteration);
    Insist(!ierr, "Error getting iteration number.");

  }

  // Diagnostic output
  if (d_print_level > 0)
  {
    printf(" MGSolverGMRES Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, norm_resid, d_sweeper->number_sweeps());
  }
  if (norm_resid > d_tolerance)
  {
    detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
      "    MGSolverGMRES did not converge.");
  }


} // end solve

// Implementation

template <class D>
inline void MGSolverGMRES<D>::build_rhs(State::vec_moments_type &B)
{
  using detran_utilities::norm_residual;

  for (int g = d_upscatter_cutoff; g < d_number_groups; g++)
  {
    int g_index = g - d_upscatter_cutoff;

    // Setup the sweeper.
    d_sweeper->setup_group(g);

    // Clear the group boundary and set with any fixed conditions.
    d_boundary->clear(g);
    d_boundary->set(g);

    // Build the fixed source.  For fixed source problems, this
    // includes external sources and any in-scatter from groups
    // solved by Gauss-Seidel.  For eigenproblems, fission is
    // also included.
    d_sweepsource->reset();
    d_sweepsource->build_fixed_with_downscatter(g, d_upscatter_cutoff);

    // If no reflective, one sweep and out.
    if (!d_boundary->has_reflective())
    {
      d_sweeper->sweep(B[g_index]);
    }
    // Otherwise, solve the reflecting condition problem.
    else
    {
      // Create dummy moment container for convergence.
      State::moments_type B_old(B[g_index].size(), 0.0);
      // Set sweeper to update boundaries on the fly.
      d_sweeper->set_update_boundary(true);
      int iteration;
      double error = 0;
      for (iteration = 0; iteration < 1000; iteration++)
      {
        // Update boundary.  This updates boundaries due to reflection, etc.
        d_boundary->update(g);
        // Swap old and new.
        std::swap(B[g_index], B_old);
        // Sweep.
        d_sweeper->sweep(B[g_index]);
        // Flux residual using L-infinity.
        error = norm_residual(B_old, B[g_index], "Linf");
        if (error < d_tolerance)
        {
          d_reflective_solve_iterations = iteration + 1;
          break;
        }
      }
      d_sweeper->set_update_boundary(false);
    }

  } // group loop

  // Fill the PETSc Vec with B. Remember to replace it.
  double *B_a;
  VecSet(d_B, 0.0);
  VecGetArray(d_B, &B_a);
  for (int g = 0; g < d_upscatter_size; g++)
  {
    int offset = g * (d_moments_size_group + d_boundary_size_group);
    for (int i = 0; i < d_moments_size_group; i++)
    {
      B_a[i + offset] = B[g][i];
    }
  }
  VecRestoreArray(d_B, &B_a);

}

template <class D>
inline PetscErrorCode MGSolverGMRES<D>::apply_MGTO(Mat A, Vec X, Vec Y)
{

  using std::cout;
  using std::endl;
  PetscErrorCode ierr;

  // Get the arrays from the Krylov vectors.
  double *X_a;
  ierr = VecGetArray(X, &X_a); CHKERRQ(ierr);
  double *Y_a;
  ierr = VecGetArray(Y, &Y_a); CHKERRQ(ierr);

//  // Assing array to phi_vac.
//  State::vec_moments_type
//    phi_original(d_upscatter_size,
//                 State::moments_type(d_moments_size_group, 0.0));
//  State::vec_moments_type
//    phi_update(d_upscatter_size,
//               State::moments_type(d_moments_size_group, 0.0));

  State::vec_moments_type phi_original(d_state->all_phi());
  State::vec_moments_type phi_update(d_state->all_phi());

  for (int g = d_upscatter_cutoff; g < d_number_groups; g++)
  {
    int offset = (g - d_upscatter_cutoff) *
                 (d_moments_size_group + d_boundary_size_group);
    for (int i = 0; i < d_moments_size_group; i++)
    {
      phi_original[g][i] = X_a[i + offset];
      phi_update[g][i]   = X_a[i + offset];
    }
  }

  // Sweep each group
  for (int g = d_upscatter_cutoff; g < d_number_groups; g++)
  {

    int g_index = g - d_upscatter_cutoff;
    int g_size  = d_moments_size_group + d_boundary_size_group;
    int m_offset = g_index * g_size;
    int b_offset = m_offset + d_moments_size_group;

    // Reset the source and place the original outgoing boundary flux.
    d_boundary->clear(g);

    if (d_boundary->has_reflective())
    {
      // Set the incident boundary flux.
      d_boundary->set_incident(g, X_a + b_offset);
    }

    // Reset the source to zero.
    d_sweepsource->reset();
    d_sweepsource->build_total_scatter(g, d_upscatter_cutoff, phi_original);

    // Set the sweeper and sweep.
    d_sweeper->setup_group(g);
    d_sweeper->sweep(phi_update[g]);

    // Update outgoing vector
    {

      // Assign the moment values.
      for (int i = 0; i < d_moments_size_group; i++)
      {
        Y_a[i + m_offset] = phi_original[g][i] - phi_update[g][i];
      }

      if (d_boundary->has_reflective())
      {
        // Update the boundary and fetch.
        d_boundary->update(g);
        vec_dbl psi_update(d_boundary_size_group, 0.0);
        d_boundary->get_incident(g, &psi_update[0]);

        // Add the boundary values.
        for (int a = 0; a < d_boundary_size_group; a++)
        {
          Y_a[a + b_offset] = X_a[a + b_offset] - psi_update[a];
        }
      }

    } // end update

  }

  // Restore the arrays.
  ierr = VecRestoreArray(X, &X_a); CHKERRQ(ierr);
  ierr = VecRestoreArray(Y, &Y_a); CHKERRQ(ierr);

  // No error.
  return 0;

}

} // end namespace detran

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

inline PetscErrorCode apply_MGTO_1D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as MGSolverGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::MGSolverGMRES<detran::_1D> *outer =
    (detran::MGSolverGMRES<detran::_1D>*) ctx;
  // Call the actual apply operator.
  return outer->apply_MGTO(A, x, y);
}

inline PetscErrorCode apply_MGTO_2D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as MGSolverGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::MGSolverGMRES<detran::_2D> *outer =
    (detran::MGSolverGMRES<detran::_2D>*) ctx;
  // Call the actual apply operator.
  return outer->apply_MGTO(A, x, y);
}

inline PetscErrorCode apply_MGTO_3D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as MGSolverGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::MGSolverGMRES<detran::_3D> *outer =
    (detran::MGSolverGMRES<detran::_3D>*) ctx;
  // Call the actual apply operator.
  return outer->apply_MGTO(A, x, y);
}

#endif /* detran_MGSOLVERGMRES_I_HH_ */

//---------------------------------------------------------------------------//
//              end of MGSolverGMRES.i.hh
//---------------------------------------------------------------------------//
