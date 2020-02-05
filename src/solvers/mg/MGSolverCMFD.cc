//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSolverCMFD.cc
 *  @brief MGSolverCMFD
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "solvers/mg/MGSolverCMFD.hh"
#include "solvers/wg/WGSolverSI.hh"
#include "utilities/MathUtilities.hh"
#include "transport/SweepSource.hh"
#include "callow/solver/LinearSolverCreator.hh"

#include <iostream>

#define COUT(c) std::cout << c << std::endl;

namespace detran
{

using detran_utilities::norm;
using detran_utilities::norm_residual;
using detran_utilities::range;
using detran_utilities::vec_scale;

//----------------------------------------------------------------------------//
template <class D>
MGSolverCMFD<D>::MGSolverCMFD(SP_state                  state,
                  SP_material               material,
                  SP_boundary               boundary,
                  const vec_externalsource &q_e,
                  SP_fissionsource          q_f,
                  bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
  , d_lower(d_adjoint ? d_number_groups - 1 : 0)
  , d_upper(d_adjoint ? -1 : d_material->number_groups())
  , d_diff_coef_weight(Homogenize::PHI_D)
  , d_omega(1.0)
{
  // Force use of source iteration
  d_wg_solver = new WGSolverSI<D>(d_state, d_material, d_quadrature,
                                  d_boundary, d_externalsources,
                                  d_fissionsource, d_multiply);

  // Create coarse mesh
  size_t level = 2;
  if (d_input->check("mgpc_coarse_mesh_level"))
  {
    level = d_input->template get<int>("mgpc_coarse_mesh_level");
    Assert(level > 0);
  }
  d_coarsener = new CoarseMesh(d_mesh, level);
  d_coarse_mesh = d_coarsener->get_coarse_mesh();

  // Create the current tally
  d_tally = new CurrentTally<D>(d_coarsener, d_quadrature, d_number_groups);
  d_wg_solver->get_sweeper()->set_tally(d_tally);

  // Check solver db
  if (d_input->check("outer_pc_db"))
    d_solver_db = d_input->template get<SP_input>("outer_pc_db");

  if (d_input->check("cmfd_diff_coef_weight"))
    d_diff_coef_weight = d_input->template get<int>("cmfd_diff_coef_weight");
  Assert(d_diff_coef_weight < Homogenize::END_OPTIONS_DC);


  if (d_input->check("cmfd_relaxation"))
    d_omega = d_input->template get<double>("cmfd_relaxation");

  // Create the linear solver for this group.
  d_solver = callow::LinearSolverCreator::Create(d_solver_db);

}

//----------------------------------------------------------------------------//
template <class D>
void MGSolverCMFD<D>::solve(const double keff)
{
  // Save current group flux
  State::group_moments_type phi_old = d_state->all_phi();

  // Norm of the group-wise residuals and the total residual norm
  vec_dbl nres(d_number_groups, 0.0);
  double nres_tot = 0;

  // Set the scaling factor for multiplying problems
  if (d_multiply) d_fissionsource->setup_outer(1.0 / keff);

  vec_size_t groups = range<size_t>(d_lower, d_upper);
  vec_size_t::iterator g_it;
  size_t iteration;
  for (iteration = 0; iteration <= d_maximum_iterations; ++iteration)
  {
    detran_utilities::vec_scale(nres, 0.0);

    // Save current group flux.
    State::group_moments_type phi_old = d_state->all_phi();

    // Loop over iteration block
    energy_sweep();
    for (g_it = groups.begin(); g_it != groups.end(); ++g_it)
    {
      nres[*g_it] = norm_residual(d_state->phi(*g_it), phi_old[*g_it], "Linf");
    }
    nres_tot = norm(nres, "Linf");

    if (d_print_level > 1  && iteration % d_print_interval == 0)
    {
      printf("  CMFD Iter: %3i  Error: %12.9f \n", (int)iteration, nres_tot);
    }
    if (nres_tot < d_tolerance) break;

    // Perform CMFD update
    update(keff);

  } // end upscatter iterations

  // Diagnostic output
  if (d_print_level > 0)
  {
    printf("  CMFD Final: Number Iters: %3i  Error: %12.9f  Sweeps: %6i \n",
           (int)iteration, nres_tot, number_sweeps());
  }

}

//----------------------------------------------------------------------------//
template <class D>
void MGSolverCMFD<D>::energy_sweep()
{
  vec_size_t groups = range<size_t>(d_lower, d_upper);
  vec_size_t::iterator g_it;
  for (g_it = groups.begin(); g_it != groups.end(); ++g_it)
  {
    d_wg_solver->solve(*g_it);
  }
}

//----------------------------------------------------------------------------//
template <class D>
void MGSolverCMFD<D>::update(const double keff)
{
  // Compute the current, if applicable.
  if (d_diff_coef_weight == Homogenize::CURRENT_D)
    compute_current();

  // Homogenize the material
  Homogenize H(d_material, d_diff_coef_weight);
  SP_material cmat = H.homogenize(d_state, d_mesh, "COARSEMESH");
  const vec2_dbl &phi = H.coarse_mesh_flux();

  // Coarse mesh diffusion operator
  d_operator = new CMFDLossOperator<D>(d_input,
                                       cmat,
                                       d_coarse_mesh,
                                       d_tally,
                                       d_multiply,
                                       d_adjoint,
                                       keff);
  d_operator->construct(phi);

  // Construct source
  callow::Vector x(d_number_groups * d_coarse_mesh->number_cells(), 0.0);
  callow::Vector b(d_number_groups * d_coarse_mesh->number_cells(), 0.0);
  const vec_int &cmap = d_mesh->mesh_map("COARSEMESH");
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    d_wg_solver->get_sweepsource()->reset();
    d_wg_solver->get_sweepsource()->build_fixed(g);
    const State::moments_type &q =
      d_wg_solver->get_sweepsource()->fixed_group_source();
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t ci = cmap[i] + g * d_coarse_mesh->number_cells();
      double v_ratio = d_mesh->volume(i) / d_coarse_mesh->volume(cmap[i]);
      b[ci] += q[i] * v_ratio;
      x[ci] += d_state->phi(g)[i] * v_ratio;
    }
  }

  callow::Vector x0(x);

  // Solve the linear system
  d_solver->set_operators(d_operator, d_solver_db);
  d_solver->solve(b, x);

  // Scale the fluxes
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t ci = cmap[i] + g * d_coarse_mesh->number_cells();
      d_state->phi(g)[i] *= (d_omega * x[ci] / x0[ci] + (1-d_omega));
    }
  }

}

//----------------------------------------------------------------------------//
template <class D>
void MGSolverCMFD<D>::compute_current()
{
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      double J[3] = {0.0, 0.0, 0.0};
      for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
      {
        for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
        {
          double w = d_quadrature->weight(a);
          double psi = d_state->psi(g, o, a)[i];
          for (size_t d = 0; d < D::dimension; ++d)
          {
            J[d] += w * psi * d_quadrature->cosines(d)[a];
          }
        }
      }
      d_state->current(g)[i] = std::sqrt(J[0]*J[0] + J[1]*J[1] + J[2]*J[2]);
    }
  }
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class MGSolverCMFD<_1D>;
template class MGSolverCMFD<_2D>;
template class MGSolverCMFD<_3D>;


} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGSolverCMFD.cc
//----------------------------------------------------------------------------//
