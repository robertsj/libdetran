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

namespace detran
{

template <class D>
void PowerIteration<D>::solve()
{

  std::cout << "Starting PI." << d_max_iters << std::endl;

  // K-eigenvalue
  double keff = 1.0;
  double keff_old = 1.0;

  // Initialize the fission density
  d_fission_source->initialize();

  // Power iterations.
  int iteration;
  double error;
  for (iteration = 0; iteration < 2; iteration++)
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

    // Compute keff.
    State::moments_type fd(d_fission_source->density());
    keff_old = keff;
    double norm_fd = norm(fd);
    double norm_fd_old = norm(fd_old);
    std::cout << "  norm_fd " << norm_fd << " norm_fd_old " << norm_fd_old << std::endl;
    keff = keff_old * norm(fd, true) / norm(fd_old, true);
    THROW("lala");
    // Compute error in fission density.
    error = norm_residual(fd, fd_old, true);
    std::cout << "  PI Iter: " << iteration << " Error: " << error <<  " keff: " << keff << std::endl;
    //if (error < d_tolerance) break;

  } // eigensolver loop

  if (error > d_tolerance) warning(SOLVER_CONVERGENCE,
                                   "PowerIteration did not converge.");

  std::cout << "PI done." << std::endl;
}

// Explicit instantiations
template class PowerIteration<_1D>;
template class PowerIteration<_2D>;
template class PowerIteration<_3D>;

} // end namespace detran

#endif /* POWERITERATION_I_HH_ */
