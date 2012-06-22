//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussSeidel.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  GaussSeidel inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef GAUSSSEIDEL_I_HH_
#define GAUSSSEIDEL_I_HH_

// Utilities
#include "Warning.hh"

// System
#include <algorithm>
#include <cstdio>

namespace detran
{

template <class D>
void GaussSeidel<D>::solve()
{

  if (d_print_out > 0) std::cout << "  Starting GS." << std::endl;

  // Initial downscatter.
  for (int g = 0; g < d_number_groups; g++)
  {
    d_inner_solver->solve(g);
  }

  // Do upscatter iterations if required.  Skip these if the max_iters = 0.
  if (!d_downscatter and d_max_iters > 0 and d_number_groups > 1)
  {
    // Upscatter iterations.
    int iteration;
    double error;
    for (iteration = 1; iteration <= d_max_iters; iteration++)
    {
      // Reset the error.
      error = 0.0;

      // Save current group flux.
      State::group_moments_type phi_old = d_state->all_phi();

      // Loop over just those groups into which upscatter occurs.  I have
      // no idea whether symmetric GS (i.e. loop up and down) is better.
      for (int g = d_material->upscatter_cutoff(); g < d_number_groups; g++)
      {
        d_inner_solver->solve(g);
        // Constructing the L-inf norm piecewise.
        error = std::max(error,
                         norm_residual(d_state->phi(g), phi_old[g], "Linf"));
      }

      if (d_print_out > 1  and iteration % d_print_interval == 0)
      {
        printf("  GS Iter: %3i  Error: %12.9f \n", iteration, error);
      }
      if (error < d_tolerance) break;

    } // upscatter loop

    if (d_print_out > 0)
    {
      printf("  GS Final: Number Iter: %3i  Error: %12.9f \n",
             iteration, error);
    }

    if (error > d_tolerance)
      warning(SOLVER_CONVERGENCE, "Gauss-Seidel did not converge.");

  }

 if (d_print_out > 0) std::cout << "  GS done." << std::endl;

}

// Explicit instantiations

template class GaussSeidel<_1D>;
template class GaussSeidel<_2D>;
template class GaussSeidel<_3D>;

} // end namespace detran

#endif /* GAUSSSEIDEL_I_HH_ */
