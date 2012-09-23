//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussSeidelMG.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  GaussSeidelMG inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef GAUSSSEIDEL_I_HH_
#define GAUSSSEIDEL_I_HH_

#include "utilities/MathUtilities.hh"
#include "utilities/Warning.hh"
#include <algorithm>
#include <cstdio>

namespace detran
{

template <class D>
void GaussSeidelMG<D>::solve()
{
  using detran_utilities::norm_residual;

  if (d_print_out > 0) std::cout << "  Starting GS." << std::endl;

  // Save current group flux.
  State::group_moments_type phi_old = d_state->all_phi();

  // Norm of the residual.
  double nres = 0.0;

  // Upscatter iterations.
  int iteration = 0;

  // Initial downscatter.
  for (int g = 0; g < d_number_groups; g++)
  {
    d_inner_solver->solve(g);
  }

  // Do upscatter iterations if required.  Skip these if the max_iters = 0.
  if (!d_downscatter and d_max_iters > 0 and d_number_groups > 1)
  {
    // Upscatter iterations.
    for (iteration = 1; iteration <= d_max_iters; iteration++)
    {
      // Reset the residual norm.
      nres = 0.0;

      // Save current group flux.
      State::group_moments_type phi_old = d_state->all_phi();

      // Loop over just those groups into which upscatter occurs.  I have
      // no idea whether symmetric GS (i.e. loop up and down) is better.
      for (int g = d_material->upscatter_cutoff(); g < d_number_groups; g++)
      {
        d_inner_solver->solve(g);

        // Constructing the norm piecewise.
        double nres_g = norm_residual(d_state->phi(g), phi_old[g], d_norm_type);
        if (d_norm_type == "Linf")
          nres = std::max(nres, nres_g);
        else if (d_norm_type == "L1")
          nres += nres_g;
        else
          nres += nres_g * nres_g;
      }
      if (d_norm_type == "L2")
        nres = std::sqrt(nres);

      if (d_print_out > 1  and iteration % d_print_interval == 0)
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
  if (d_print_out > 0)
  {
    int number_sweeps = d_inner_solver->get_sweeper()->number_sweeps();
    printf("  GS Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, nres, number_sweeps);
    std::cout << "  GS done!" << std::endl;
  }


}

// Explicit instantiations

template class GaussSeidelMG<_1D>;
template class GaussSeidelMG<_2D>;
template class GaussSeidelMG<_3D>;

} // end namespace detran

#endif /* GAUSSSEIDEL_I_HH_ */
