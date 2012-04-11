//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PowerIteration.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  PowerIteration inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef POWERITERATION_I_HH_
#define POWERITERATION_I_HH_

// Utilities
#include "MathUtilities.hh"
#include "Warning.hh"

// System
#include <algorithm>
#include <cstdio>

namespace detran
{

template <class D>
void PowerIteration<D>::solve()
{

  std::cout << "Starting PI." << std::endl;

  // K-eigenvalue
  double keff = 1.0;
  double keff_old = 1.0;

  // Initialize the fission density
  d_fission_source->initialize();

  // Power iterations.
  int iteration;
  double error;
  for (iteration = 1; iteration <= d_max_iters; iteration++)
  {
    // Reset the error.
    error = 0.0;

    // Save current density.
    State::moments_type fd_old(d_fission_source->density());

    // Setup outer iteration.
    d_fission_source->setup_outer(1/keff);

    // Solve the multigroup equations.
    d_mg_solver->solve();

    // Update density.
    d_fission_source->update();

    // Compute keff.  Here, using L1.  Could implement
    // volume-integrated fission rate if desired.
    State::moments_type fd(d_fission_source->density());
    keff_old = keff;
    keff = keff_old * norm(fd, "L1") / norm(fd_old, "L1");

    // Compute error in fission density.
    error = norm_residual(fd, fd_old, "L1");

    if (d_print_out > 1 and iteration % d_print_interval == 0)
    {
      printf("PI Iter: %3i  Error: %12.9f  keff: %12.9f \n",
             iteration, error, keff);
    }
    if (error < d_tolerance) break;

  } // eigensolver loop

  if (d_print_out > 0)
  {
    printf("*********************************************************************\n");
    printf(" PI Final: Number Iters: %3i  Error: %12.9f keff: %12.9f \n",
           iteration, error, keff);
    printf("*********************************************************************\n");
  }

  if (error > d_tolerance) warning(SOLVER_CONVERGENCE,
                                   "PowerIteration did not converge.");
  d_state->set_eigenvalue(keff);
  std::cout << "PI done." << std::endl;
}

// Explicit instantiations
template class PowerIteration<_1D>;
template class PowerIteration<_2D>;
template class PowerIteration<_3D>;

} // end namespace detran

#endif /* POWERITERATION_I_HH_ */
