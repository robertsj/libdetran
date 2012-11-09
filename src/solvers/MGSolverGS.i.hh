//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGSolverGS.i.hh
 *  @author robertsj
 *  @date   Apr 10, 2012
 *  @brief  MGSolverGS inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGSolverGS_I_HH_
#define detran_MGSolverGS_I_HH_

#include "ioutils/StdOutUtils.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/Warning.hh"
#include <algorithm>
#include <cstdio>

namespace detran
{

template <class D>
void MGSolverGS<D>::solve(const double keff)
{
  using detran_utilities::norm_residual;

  if (d_print_level > 0) std::cout << "  Starting GS." << std::endl;

  // Save current group flux.
  State::group_moments_type phi_old = d_state->all_phi();

  // Norm of the residual.
  double nres = 0.0;

  // Upscatter iterations.
  int iteration = 0;

  // Set the scaling factor for multiplying problems
  if (d_multiply) d_fissionsource->setup_outer(1.0/keff);

  // Initial downscatter.
  for (int g = 0; g < d_number_groups; g++)
  {
    // Recompute the fission density
    if (d_multiply) d_fissionsource->update();
    // Solve the within-group equation
    d_wg_solver->solve(g);
  }

  // Decide whether to iterate or not
  bool iterate = false;
  if ((!d_downscatter and d_maximum_iterations > 0 and d_number_groups > 1)
      or d_multiply)
  {
    iterate = true;
  }

  // Do upscatter iterations if required.  Skip these if the max_iters = 0.
  if (iterate)
  {
    // Iterations
    for (iteration = 1; iteration <= d_maximum_iterations; iteration++)
    {

      // Iteration residual norm for all fluxes and group flux
      double nres = 0.0;
      double nres_g = 0.0;

      // Save current group flux.
      State::group_moments_type phi_old = d_state->all_phi();

      // Set group iteration lower bound
      int g_lower = d_material->upscatter_cutoff();
      if (d_multiply) g_lower = 0;

      // Loop over required groups
      for (int g = g_lower; g < d_number_groups; ++g)
      {
        // Recompute the fission density
        if (d_multiply) d_fissionsource->update();
        // Solve the within-group equation
        d_wg_solver->solve(g);
        // Constructing the norm piecewise.
        nres_g = norm_residual(d_state->phi(g), phi_old[g], d_norm_type);
        if (d_norm_type == "Linf")
          nres = std::max(nres, nres_g);
        else if (d_norm_type == "L1")
          nres += nres_g;
        else
          nres += nres_g * nres_g;
      }
      if (d_norm_type == "L2")
        nres = std::sqrt(nres);

      if (d_print_level > 1  and iteration % d_print_interval == 0)
      {
        printf("  GS Iter: %3i  Error: %12.9f \n", iteration, nres);
      }
      if (nres < d_tolerance) break;

    } // end upscatter iterations

    if (nres > d_tolerance)
    {
      detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
        "Gauss-Seidel upscatter did not converge.");
    }

  } // end upscatter block

  // Diagnostic output
  if (d_print_level > 0)
  {
    int number_sweeps = d_wg_solver->get_sweeper()->number_sweeps();
    printf("  GS Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, nres, number_sweeps);
    std::cout << "  GS done!" << std::endl;
  }

}

} // end namespace detran

#endif /* detran_MGSolverGS_I_HH_ */
