//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenPI.i.hh
 *  @brief EigenPI inline member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_EIGENPI_I_HH_
#define detran_EIGENPI_I_HH_

#include "utilities/MathUtilities.hh"
#include "utilities/Warning.hh"
#include <algorithm>
#include <cstdio>

namespace detran
{

template <class D>
void EigenPI<D>::solve()
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
  d_fissionsource->initialize();
//
//  d_state->clear();
//  for (int g = 0; g < d_number_groups; ++g)
//    for (int i = 0; i < d_state->get_mesh()->number_cells(); ++i)
//      d_state->phi(g)[i] = 1.0;
//  d_fissionsource->update();

  // Power iterations.
  int iteration;
  double error;

  for (iteration = 1; iteration <= d_maximum_iterations; iteration++)
  {
    // Reset the error.
    error = 0.0;

    // Save current density.
    State::moments_type fd_old(d_fissionsource->density());

    // Setup outer iteration.  This precomputes the group sources.
    d_fissionsource->setup_outer(1/keff);

    // Solve the multigroup equations.
    d_mg_solver->solve();

    // Update density.
    d_fissionsource->update();

    // Compute keff.  Here, using L1.  Could implement
    // volume-integrated fission rate if desired.
    State::moments_type fd(d_fissionsource->density());

    // Overrelaxation
    if (d_omega != 1.0)
    {
      for (int i = 0; i < fd.size(); ++i)
        fd[i] = d_omega * fd[i] + (1.0 - d_omega) * fd_old[i];
    }
    keff_2 = keff_1;
    keff_1 = keff;
    keff   = keff_1 * norm(fd, "L1") / norm(fd_old, "L1");

    // Compute error in fission density.
    error = norm_residual(fd, fd_old, "L1");
    if (d_print_level > 1 && iteration % d_print_interval == 0)
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
    if (error < d_tolerance) break;

  } // eigensolver loop

  if (d_print_level > 0)
  {
    printf("*********************************************************************\n");
    printf(" PI Final: Number Iters: %3i  Error: %12.9e keff: %12.9f \n",
           iteration, error, keff);
    printf("*********************************************************************\n");
  }

  if (error > d_tolerance)
  {
    detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
      "EigenPI did not converge.");
  }
  d_state->set_eigenvalue(keff);
  std::cout << "PI done." << std::endl;
}

} // end namespace detran

#endif /* detran_EIGENPI_I_HH_ */

//----------------------------------------------------------------------------//
//              end of file EigenPI.i.hh
//----------------------------------------------------------------------------//

