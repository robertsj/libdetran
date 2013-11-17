//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverGS.i.hh
 *  @brief MGSolverGS inline member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGSolverGS_I_HH_
#define detran_MGSolverGS_I_HH_

#include "ioutils/StdOutUtils.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/Warning.hh"
#include <algorithm>
#include <cstdio>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
void MGSolverGS<D>::solve(const double keff)
{
  using detran_utilities::norm;
  using detran_utilities::norm_residual;
  using detran_utilities::range;
  using detran_utilities::vec_scale;

  // Save current group flux
  State::group_moments_type phi_old = d_state->all_phi();

  // Norm of the group-wise residuals and the total residual norm
  vec_dbl nres(d_number_groups, 0.0);
  double nres_tot = 0;

  // Set the scaling factor for multiplying problems
  if (d_multiply) d_fissionsource->setup_outer(1.0 / keff);

  // Initial sweep through all groups.
  vec_size_t groups = range<size_t>(d_lower, d_upper);
  vec_size_t::iterator g_it = groups.begin();
  for (; g_it != groups.end(); ++g_it)
  {
    d_wg_solver->solve(*g_it);
  }

  // Perform upscatter iterations.
  size_t iteration = 0;
  if (d_iterate)
  {
    // Update group bounds for the possibly truncated iteration block.
    groups = range<size_t>(d_lower_upscatter, d_upper);

    for (iteration = 1; iteration <= d_maximum_iterations; ++iteration)
    {
      detran_utilities::vec_scale(nres, 0.0);

      // Save current group flux.
      State::group_moments_type phi_old = d_state->all_phi();

      // Loop over iteration block
      for (g_it = groups.begin(); g_it != groups.end(); ++g_it)
      {
        size_t g = *g_it;
        //if (d_multiply) d_fissionsource->update();
        d_wg_solver->solve(g);
        nres[g] = norm_residual(d_state->phi(g), phi_old[g], d_norm_type);
      }
      nres_tot = norm(nres, d_norm_type);

      if (d_print_level > 1  && iteration % d_print_interval == 0)
      {
        printf("  GS Iter: %3i  Error: %12.9f \n", (int)iteration, nres_tot);
      }
      if (nres_tot < d_tolerance) break;

    } // end upscatter iterations

    if (nres_tot > d_tolerance)
    {
      detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
        "Gauss-Seidel upscatter did not converge.");
    }

  } // end upscatter block

  // Diagnostic output
  if (d_print_level > 0)
  {
    printf("  GS Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           (int)iteration, nres_tot, number_sweeps());
  }

}

//----------------------------------------------------------------------------//
template <class D>
void MGSolverGS<D>::sweep()
{
  using detran_utilities::range;
  vec_size_t groups = range<size_t>(d_lower, d_upper);
  vec_size_t::iterator g_it = groups.begin();
  for (; g_it != groups.end(); ++g_it)
  {
    d_wg_solver->solve(*g_it);
  }
}

} // end namespace detran

#endif /* detran_MGSolverGS_I_HH_ */

//----------------------------------------------------------------------------//
//              end of MGSolverGS.i.hh
//----------------------------------------------------------------------------//
