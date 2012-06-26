//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SourceIteration.i.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  SourceIteration inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef SOURCEITERATION_I_HH_
#define SOURCEITERATION_I_HH_

// System
#include <algorithm>
#include <cstdio>
#include <iostream>

namespace detran
{

template <class D>
void SourceIteration<D>::solve(int g)
{
  using std::cout;
  using std::endl;

  if (d_print_out > 0) std::cout << "    Starting SI." << std::endl;

  // Setup boundary conditions.  This sets any conditions fixed for the solve.
  d_boundary->set(g);

  // Set the equations.
  d_sweeper->setup_group(g);

  // Build the fixed source.  This includes external, fission, and in-scatter.
  d_sweepsource->build_fixed_with_scatter(g);

  // Get flux and allocate last flux.
  moments_type phi(d_state->phi(g));
  moments_type phi_old(phi.size(), 0.0);

  // Construct within group.
  d_sweepsource->build_within_group_scatter(g, phi);

  // Iterate.
  double error = 1.0;
  int iteration;
  for (iteration = 1; iteration <= d_max_iters; iteration++)
  {

    // Update boundary.  This updates boundaries due to reflection, etc.
    //d_boundary->update(g); // sweeper should update this.

    // Swap old and new.
    std::swap(phi, phi_old);

    // Sweep.
    d_sweeper->sweep(phi);

    // Flux residual using L-infinity.
    error = norm_residual(phi_old, phi, "Linf");

    if (d_print_out > 1 and iteration % d_print_interval == 0)
    {
      printf("    SI Iter: %3i  Error: %12.9f \n", iteration, error);
    }
    if (error < d_tolerance) break;

    // Construct within group
    d_sweepsource->build_within_group_scatter(g, phi);

  } // end iterations

  if (d_print_out > 0)
  {
    printf("    SI Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, error, d_sweeper->number_sweeps());
  }

  if (error > d_tolerance)
    warning(SOLVER_CONVERGENCE, "    SourceIteration did not converge.");

  // Update the state with the new flux.
  d_state->phi(g) = phi;

}

template class SourceIteration<_1D>;
template class SourceIteration<_2D>;
template class SourceIteration<_3D>;

} // namespace detran

#endif /* SOURCEITERATION_I_HH_ */

//---------------------------------------------------------------------------//
//              end of SourceIteration.i.hh
//---------------------------------------------------------------------------//
