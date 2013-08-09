//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverGMRES.i.hh
 *  @brief MGSolverGMRES inline member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGSOLVERGMRES_I_HH_
#define detran_MGSOLVERGMRES_I_HH_

#include "utilities/MathUtilities.hh"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cstring>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
inline void MGSolverGMRES<D>::solve(const double keff)
{
  using detran_utilities::range;

  double norm_resid = 0.0;
  int iteration = 0;

  // Set the scaling factor for multiplying problems
  if (d_multiply) d_fissionsource->set_scale(1.0 / keff);

  // Build the preconditioner
  if (d_pc) d_pc->build(keff, d_state);

  // Debug printing of operators
  bool flag = false;
  if (d_input->check("print_transport_operator"))
  {
    flag = d_input->template get<int>("print_transport_operator") != 0;
    if (flag != 0) d_operator->compute_explicit("tran.out");
  }
  if (d_input->check("print_preconditioner_operator"))
  {
    flag = d_input->template get<int>("print_preconditioner_operator") != 0;
    if (flag && d_pc) d_pc->display("pc.out");
  }

  //--------------------------------------------------------------------------//
  // GAUSS-SEIDEL BLOCK
  //--------------------------------------------------------------------------//

  groups_t groups  = range<size_t>(d_lower, d_krylov_group_cutoff);
  groups_iter g_it = groups.begin();
  for (; g_it != groups.end(); ++g_it)
  {
    d_wg_solver->solve(*g_it);
  }

  //--------------------------------------------------------------------------//
  // KRYLOV BLOCK
  //--------------------------------------------------------------------------//

  if (d_number_active_groups)
  {
    groups = range<size_t>(d_krylov_group_cutoff, d_upper);

    //------------------------------------------------------------------------//
    // BUILD RIGHT HAND SIDE
    //------------------------------------------------------------------------//

    State::vec_moments_type B(d_number_active_groups,
                              State::moments_type(d_moments_size_group, 0.0));
    build_rhs(B);

    //------------------------------------------------------------------------//
    // SOLVE MULTIGROUP TRANSPORT EQUATION
    //------------------------------------------------------------------------//

    // Set the uncollided flux as the initial guess.  It may be of value to
    // have an optional initial guess based on the current flux.
    d_x->copy(d_b);
    d_solver->solve(*d_b, *d_x);

    //------------------------------------------------------------------------//
    // POST-PROCESS
    //------------------------------------------------------------------------//

    // Copy the flux solution into state.
    for (g_it = groups.begin(); g_it != groups.end(); ++g_it)
    {
      int g_index = d_adjoint ? *g_it : *g_it - d_krylov_group_cutoff;
      int offset  = g_index * d_moments_size_group;
      memcpy(&d_state->phi(*g_it)[0],
             &((*d_x)[offset]),
             d_moments_size_group*sizeof(double));
    }

    // Set the incident boundary flux if applicable
    if (d_boundary->has_reflective())
    {
      for (g_it = groups.begin(); g_it != groups.end(); ++g_it)
      {
        int g_index = d_adjoint ? *g_it : *g_it - d_krylov_group_cutoff;
        int offset  = g_index * d_boundary_size_group +
                      d_moments_size_group * d_number_active_groups;
        d_boundary->psi(*g_it, &(*d_x)[offset],
                        BoundaryBase<D>::IN, BoundaryBase<D>::SET, true);
      }
    }

    // Iterate to update the angular fluxes.  This includes cell and
    // boundary fluxes.  Remember, the angular fluxes are only a secondary
    // output of the solvers and are not part of the Krylov solution vector.
    if (d_update_angular_flux)
    {
      for (g_it = groups.begin(); g_it != groups.end(); ++g_it)
      {
        unsigned int g = *g_it;
        d_boundary->set(g);
        State::moments_type phi_g = d_state->phi(g);
        d_sweeper->setup_group(g);
        d_sweepsource->reset();
        d_sweepsource->build_fixed_with_scatter(g);
        d_sweepsource->build_within_group_scatter(g, phi_g);
        solve_wg_reflection(g, phi_g);
      }
    }

  } // end krylov block

  // Diagnostic output
  if (d_print_level > 0)
  {
    printf(" MGSolverGMRES Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, norm_resid, number_sweeps());
  }
  if (norm_resid > d_tolerance)
  {
    detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
      "    MGSolverGMRES did not converge.");
  }


} // end solve

//----------------------------------------------------------------------------//
template <class D>
inline void MGSolverGMRES<D>::build_rhs(State::vec_moments_type &B)
{
  // Perform a group sweep
  group_sweep(B);

  // Fill the source vector
  for (int g = 0; g < d_number_active_groups; g++)
  {
    int offset = g * d_moments_size_group;
    for (int i = 0; i < d_moments_size_group; i++)
      (*d_b)[i + offset] = B[g][i];
  }
}

//----------------------------------------------------------------------------//
template <class D>
inline void MGSolverGMRES<D>::
group_sweep(State::vec_moments_type &phi)
{
  using detran_utilities::range;
  groups_t groups  = range<size_t>(d_krylov_group_cutoff, d_upper);
  groups_iter g_it = groups.begin();
  for (; g_it != groups.end(); ++g_it)
  {
    int g = *g_it;
    int g_index = d_adjoint ? g : g - d_krylov_group_cutoff;

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
    d_sweepsource->build_fixed_with_downscatter(g, d_krylov_group_cutoff);

    // If no reflective, one sweep and out. Otherwise, iterate.
    if (!d_boundary->has_reflective())
      d_sweeper->sweep(phi[g_index]);
    else
      solve_wg_reflection(g, phi[g_index]);
  }
}

//----------------------------------------------------------------------------//
template <class D>
inline void MGSolverGMRES<D>::
solve_wg_reflection(size_t g, State::moments_type &phi_g)
{
  // Create dummy moment container for convergence.
  State::moments_type phi_g_old(phi_g.size(), 0.0);

  // Set sweeper to update boundaries on the fly.
  d_sweeper->set_update_boundary(true);

  // Iterate on the reflective conditions
  int iteration;
  double error = 0;
  for (iteration = 0; iteration < 2000; iteration++)
  {
    // Update boundary.  This updates boundaries due to reflection.
    d_boundary->update(g);

    // Swap old and new.
    std::swap(phi_g, phi_g_old);

    // Sweep
    d_sweeper->sweep(phi_g);

    // Compute residual and check convergence
    error = detran_utilities::norm_residual(phi_g_old, phi_g, "Linf");

    if (error < d_tolerance) break;
  }
  d_reflective_solve_iterations += iteration;

  // Switch sweeper back to not updating boundaries
  d_sweeper->set_update_boundary(false);
}

} // end namespace detran

#endif /* detran_MGSOLVERGMRES_I_HH_ */

//----------------------------------------------------------------------------//
//              end of MGSolverGMRES.i.hh
//----------------------------------------------------------------------------//
