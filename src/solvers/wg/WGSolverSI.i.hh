//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   WGSolverSI.i.hh
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  WGSolverSI inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_WGSOLVERSI_I_HH
#define detran_WGSOLVERSI_I_HH

// System
#include <algorithm>
#include <cstdio>
#include <iostream>

namespace detran
{

template <class D>
void WGSolverSI<D>::solve(const size_t g)
{
  using std::cout;
  using std::endl;
  using detran_utilities::norm_residual;

  if (d_print_level > 0) std::cout << "    Starting SI." << std::endl;
  if (d_print_level > 0) std::cout << "      group " << g << std::endl;

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
  for (iteration = 1; iteration <= d_maximum_iterations; iteration++)
  {

    // Update boundary.  This updates boundaries due to reflection, etc.
    //d_boundary->update(g); // sweeper should update this.

    // Swap old and new.
    std::swap(phi, phi_old);

    // Sweep.
    d_sweeper->sweep(phi);

    // Flux residual using L-infinity.
    error = norm_residual(phi_old, phi, "Linf");

    if (d_print_level > 1 and iteration % d_print_interval == 0)
    {
      printf("    SI Iter: %3i  Error: %12.9f \n", iteration, error);
    }
    if (error < d_tolerance) break;

    // Construct within group
    d_sweepsource->build_within_group_scatter(g, phi);

  } // end iterations

  if (d_print_level > 0)
  {
    printf("    SI Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           iteration, error, d_sweeper->number_sweeps());
  }

//  if (error > d_tolerance)
//  {
//    detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
//      "    WGSolverSI did not converge.");
//  }

  // Update the state with the new flux.
  d_state->phi(g) = phi;

}

} // namespace detran

#endif /* detran_WGSOLVERSI_I_HH */

//---------------------------------------------------------------------------//
//              end of WGSolverSI.i.hh
//---------------------------------------------------------------------------//
