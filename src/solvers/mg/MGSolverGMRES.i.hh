//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGSolverGMRES.i.hh
 *  @author robertsj
 *  @date   Jun 19, 2012
 *  @brief  MGSolverGMRES inline member definitions.
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

//---------------------------------------------------------------------------//
template <class D>
inline void MGSolverGMRES<D>::solve(const double keff)
{
  using std::cout;
  using std::endl;


  double norm_resid = 0.0;
  int iteration = 0;

  //-------------------------------------------------------------------------//
  // GAUSS-SEIDEL BLOCK
  //-------------------------------------------------------------------------//

  for (int g = 0; g < d_krylov_group_cutoff; g++)
  {
    d_wg_solver->solve(g);
  }

  //-------------------------------------------------------------------------//
  // KRYLOV BLOCK
  //-------------------------------------------------------------------------//

  if (d_number_active_groups)
  {

    //-----------------------------------------------------------------------//
    // BUILD RIGHT HAND SIDE
    //-----------------------------------------------------------------------//

    State::vec_moments_type B(d_number_active_groups,
                              State::moments_type(d_moments_size_group, 0.0));
    build_rhs(B);

    //-----------------------------------------------------------------------//
    // SOLVE MULTIGROUP TRANSPORT EQUATION
    //-----------------------------------------------------------------------//

    d_x->copy(d_b);
    d_solver->solve(*d_b, *d_x);

    //-----------------------------------------------------------------------//
    // POST-PROCESS
    //-----------------------------------------------------------------------//

    // Copy the flux solution into state.
    for (int g = d_krylov_group_cutoff; g < d_number_groups; g++)
    {
      int offset = (g - d_krylov_group_cutoff) * d_moments_size_group;
      memcpy(&d_state->phi(g)[0],
             &((*d_x)[offset]),
             d_moments_size_group*sizeof(double));
    }

    // Set the incident boundary flux if applicable
    if (d_boundary->has_reflective())
    {
      for (int g = d_krylov_group_cutoff; g < d_number_groups; g++)
      {
        int offset = (g - d_krylov_group_cutoff) * d_boundary_size_group +
                     d_moments_size_group * d_number_active_groups;
        d_boundary->psi(g, &(*d_x)[offset],
                        BoundaryBase<D>::IN, BoundaryBase<D>::SET, true);
      }
    }

    // Iterate to pick up outgoing boundary fluxes.  Only do this if
    // there is at least one vacuum boundary.  Otherwise, the
    // boundary fluxes are never going to be needed.
    if (d_update_boundary_flux)
    {
      //d_boundary->clear_bc();
      for (int g = d_krylov_group_cutoff; g < d_number_groups; g++)
      {
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
           iteration, norm_resid, d_sweeper->number_sweeps());
  }
  if (norm_resid > d_tolerance)
  {
    detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
      "    MGSolverGMRES did not converge.");
  }


} // end solve

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
template <class D>
inline void MGSolverGMRES<D>::
group_sweep(State::vec_moments_type &phi)
{
  for (int g = d_krylov_group_cutoff; g < d_number_groups; g++)
  {
    int g_index = g - d_krylov_group_cutoff;

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

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
//              end of MGSolverGMRES.i.hh
//---------------------------------------------------------------------------//
