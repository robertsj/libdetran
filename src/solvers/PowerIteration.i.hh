//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIteration.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  PowerIteration inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef POWERITERATION_I_HH_
#define POWERITERATION_I_HH_

#include "utilities/MathUtilities.hh"
#include "utilities/Warning.hh"
#include <algorithm>
#include <cstdio>

namespace detran
{

template <class D>
void PowerIteration<D>::solve()
{
  using detran_utilities::norm_residual;
  using detran_utilities::norm;

  std::cout << "Starting PI." << std::endl;

  // New k-eigenvalue
  double keff = 1.0;
  // k-eigenvalue from 1 time ago
  double keff_1 = 1.0;
  // k-eigenvalue from 2 times ago
  double keff_2 = 1.0;

  // Initialize the fission density
  b_fissionsource->initialize();

  // Power iterations.
  int iteration;
  double error;
  for (iteration = 1; iteration <= b_max_iters; iteration++)
  {
    // Reset the error.
    error = 0.0;

    // Save current density.
    State::moments_type fb_old(b_fissionsource->density());

    // Setup outer iteration.
    b_fissionsource->setup_outer(1/keff);

    // Solve the multigroup equations.
    b_mg_solver->solve();

    // Update density.
    b_fissionsource->update();

    // Compute keff.  Here, using L1.  Could implement
    // volume-integrated fission rate if desired.
    State::moments_type fd(b_fissionsource->density());
    keff_2 = keff_1;
    keff_1 = keff;
    keff = keff_1 * norm(fd, "L1") / norm(fb_old, "L1");

    // Compute error in fission density.
    error = norm_residual(fd, fb_old, "L1");
    if (b_print_out > 1 and iteration % b_print_interval == 0)
    {
      if (d_aitken)
      {
        double keff_aitken =
            keff_2 - (keff_1 - keff_2)*(keff_1 - keff_2)/
                     (keff - 2*keff_1 + keff_2);
        printf("PI Iter: %3i  Error: %12.9f  keff: %12.9f keffa: %12.9f \n",
               iteration, error, keff, keff_aitken);
      }
      else
      {
        printf("PI Iter: %3i  Error: %12.9f  keff: %12.9f \n",
               iteration, error, keff);
      }
    }
    if (error < b_tolerance) break;

  } // eigensolver loop

  if (b_print_out > 0)
  {
    printf("*********************************************************************\n");
    printf(" PI Final: Number Iters: %3i  Error: %12.9f keff: %12.9f \n",
           iteration, error, keff);
    printf("*********************************************************************\n");
  }

  if (error > b_tolerance)
  {
    detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
      "PowerIteration did not converge.");
  }
  b_state->set_eigenvalue(keff);
  std::cout << "PI done." << std::endl;
}

// Explicit instantiations
template class PowerIteration<_1D>;
template class PowerIteration<_2D>;
template class PowerIteration<_3D>;

} // end namespace detran

#endif /* POWERITERATION_I_HH_ */
