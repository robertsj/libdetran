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

namespace detran
{

template <class D>
void SourceIteration<D>::solve(int g)
{
  // Make some things visible.
//  using InnerIteration<D>::d_sweeper;
//  using SourceIteration<D>::d_sweepsource;
//  using SourceIteration<D>::d_max_iters;

  // Setup boundary conditions.

  // Set the equations.
  d_sweeper->setup_group(g);

  // Build the fixed source.
  d_sweepsource->build_fixed_with_scatter(g);

  // Get flux and allocate last flux.
  moments_type phi(d_state->phi(g));
  moments_type phi_old(phi.size(), 0.0);

  // Iterate.
  double residual = 0.0;
  int iter;
  for (iter = 0; iter < d_max_iters; iter++)
  {

    // Construct within group
    d_sweepsource->build_within_group_scatter(g, phi);

    // Update boundary

    // Sweep
    d_sweeper->sweep(phi);

    // Flux residual
    residual = detran_utils::norm(phi, phi_old);
    if (residual < d_tolerance)
    {
      break;
    }

    // Store flux without copying.
    std::swap(phi, phi_old);

  } // end iterations

  // did we converge?

  // Replace flux.
  d_state->phi(g) = phi;

}


} // namespace detran

#endif /* SOURCEITERATION_I_HH_ */

//---------------------------------------------------------------------------//
//              end of SourceIteration.i.hh
//---------------------------------------------------------------------------//
