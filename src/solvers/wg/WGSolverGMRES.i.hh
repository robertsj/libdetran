//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   WGSolverGMRES.i.hh
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  WGSolverGMRES inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_WGSOLVERGMRES_I_HH_
#define detran_WGSOLVERGMRES_I_HH_

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cstring>

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
inline void WGSolverGMRES<D>::solve(const size_t g)
{
  using std::cout;
  using std::endl;

  if (d_print_level > 0) std::cout << "    Starting GMRES." << std::endl;
  if (d_print_level > 0) std::cout << "      group " << g << std::endl;

  // Set the group for this solve.
  d_g = g;
  d_sweeper->setup_group(g);
  d_operator->set_group(g);
  if (d_pc) d_pc->set_group(g);

  //-------------------------------------------------------------------------//
  // BUILD RIGHT HAND SIDE
  //-------------------------------------------------------------------------//

  State::moments_type B(d_mesh->number_cells(), 0.0);
  build_rhs(B);
  double b_norm = d_b->norm(callow::L1);
  d_b->print_matlab("b.out");
  d_operator->compute_explicit("A.out");
  //d_operator->multiply(*d_b, *d_x);
  //d_x->print_matlab("x0.out");

  //-------------------------------------------------------------------------//
  // SOLVE THE TRANSPORT EQUATION
  //-------------------------------------------------------------------------//

  // Start with the uncollided flux as the initial guess
  d_x->copy(d_b);

  // Solve
  if (b_norm > 0.0) d_solver->solve(*d_b, *d_x);
  //d_x->print_matlab("x.out");

  //-------------------------------------------------------------------------//
  // POSTPROCESS
  //-------------------------------------------------------------------------//

  // Copy the flux solution into state.
  memcpy(&d_state->phi(g)[0], &(*d_x)[0],
         d_state->moments_size()*sizeof(double));

  // Set the incident boundary flux if applicable
  if (d_boundary->has_reflective())
  {
    d_boundary->psi(g, &(*d_x)[d_state->moments_size()],
                    BoundaryBase<D>::IN, BoundaryBase<D>::SET, true);
  }

  // Iterate to pick up outgoing boundary fluxes if requested.
  if (d_update_boundary_flux)
  {
    moments_type phi_g   = d_state->phi(g);
    moments_type phi_g_o = d_state->phi(g);
    d_sweeper->set_update_boundary(true);
    d_sweepsource->reset();
    d_sweepsource->build_fixed_with_scatter(d_g);
    d_sweepsource->build_within_group_scatter(d_g, phi_g);
    int iteration;
    double error = 0;
    for (iteration = 0; iteration < 1000; iteration++)
    {
      d_boundary->update(d_g);
      std::swap(phi_g, phi_g_o);
      d_sweeper->sweep(phi_g);
      error = detran_utilities::norm_residual(phi_g_o, phi_g, "Linf");
      if (error < 1e-9) break;
    }
    d_reflective_solve_iterations = iteration;
    d_sweeper->set_update_boundary(false);
  }
  //d_state->display();

  if (d_print_level > 0)
  {
    printf(" GMRES Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           d_solver->number_iterations(),
           d_solver->residual_norms()[d_solver->number_iterations()],
           d_sweeper->number_sweeps());
    if (d_solver->residual_norms()[d_solver->number_iterations()] > d_tolerance)
    {
      detran_utilities::warning(detran_utilities::SOLVER_CONVERGENCE,
        "    WGSolverGMRES did not converge.");
    }
  }

}

//---------------------------------------------------------------------------//
template <class D>
inline void WGSolverGMRES<D>::build_rhs(State::moments_type &B)
{
  using detran_utilities::norm_residual;

  d_b->set(0.0);

  // Clear the group boundary.  This is because the right hand side should
  // be based only on fixed boundaries.
  if (!d_use_initial_boundary_flux) d_boundary->clear(d_g);

  // Setup boundary conditions.  This sets any conditions fixed for the solve.
  d_boundary->set(d_g);

  // Build the fixed source.  This includes external, fission, and in-scatter.
  d_sweepsource->reset();
  d_sweepsource->build_fixed_with_scatter(d_g);

  // If no reflective, one sweep and out.
  if (!d_boundary->has_reflective())
  {
    d_sweeper->sweep(B);
  }
  // Otherwise, solve the reflecting condition problem.  This
  // needs to be redone more formally.  It performs poorly
  // for problems with lots of opposing reflection.
  else
  {
    // Create dummy moment container for convergence.
    State::moments_type B_old(B.size(), 0.0);
    // Set sweeper to update boundaries on the fly.
    d_sweeper->set_update_boundary(true);
    int iteration;
    double error = 0;
    for (iteration = 0; iteration < 1000; iteration++)
    {
      // Update boundary.  This updates boundaries due to reflection, etc.
      d_boundary->update(d_g);
      // Swap old and new.
      std::swap(B, B_old);
      // Sweep.
      d_sweeper->sweep(B);
      // Flux residual using L-infinity.
      error = norm_residual(B_old, B, "Linf");
      if (error < d_tolerance)
      {
        d_reflective_solve_iterations = iteration + 1;
        break;
      }
    }
    //std::cout << " ref its " << d_reflective_solve_iterations << std::endl;
    d_sweeper->set_update_boundary(false);
  }

  // Fill the source vector.  The RHS corresponding to the
  // boundaries is set to zero.
  for (int i = 0; i < B.size(); i++)
    (*d_b)[i] = B[i];
}

} // namespace detran

#endif /* detran_WGSOLVERGMRES_I_HH_ */

//---------------------------------------------------------------------------//
//              end of WGSolverGMRES.i.hh
//---------------------------------------------------------------------------//
