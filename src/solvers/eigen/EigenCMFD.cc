//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenCMFD.cc
 *  @brief EigenCMFD class definition
 *  @note  Copyright(C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "solvers/eigen/EigenCMFD.hh"
#include "solvers/mg/MGSolverCMFD.hh"
#include "solvers/mg/DiffusionGainOperator.hh"
#include "callow/solver/EigenSolverCreator.hh"

#include <iostream>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
EigenCMFD<D>::EigenCMFD(SP_mg_solver mg_solver)
  : Base(mg_solver)
  , d_omega(1.0)
{
  // Assert I've gotten a CMFD fixed source solver

  if (d_input->check("cmfd_relaxation"))
    d_omega = d_input->template get<double>("cmfd_relaxation");

  // Get callow solver parameter database
  if (d_input->check("eigen_solver_db"))
  {
    d_eigen_db = d_input->template get<SP_input>("eigen_solver_db");
  }
  d_eigensolver = callow::EigenSolverCreator::Create(d_eigen_db);
  Assert(d_eigensolver);

}

//----------------------------------------------------------------------------//
template <class D>
void EigenCMFD<D>::solve()
{
  using detran_utilities::norm_residual;
  using detran_utilities::norm;

  typedef MGSolverCMFD<D> Solver_T;
  Solver_T *solver = dynamic_cast<Solver_T*>(&(*d_mg_solver->solver()));
  Insist(solver, "EigenCMFD needs CMFD as the MG problem to get operator.");

  std::cout << "Starting CMFD." << std::endl;

  // New k-eigenvalue
  double keff = 1.0, keff_1 = 1.0;

  // Assume flat flux everywhere
  double val = 1.0 / d_mesh->number_cells();
  for (size_t g = 0; g < d_number_groups; ++g)
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
      d_state->phi(g)[i] = val;

  // Power iterations.
  int iteration;
  double error;

  for (iteration = 1; iteration <= d_maximum_iterations; iteration++)
  {
    // Reset the error.
    error = 0.0;

    d_fissionsource->update();
    d_fissionsource->setup_outer(1.0 / keff);
    State::moments_type fd_old(d_fissionsource->density());

    // Do one energy sweep
    solver->energy_sweep();

    // CMFD update
    double keff_c = cmfd_update();
    //keff = keff_c;

    // Update fission source.
    d_fissionsource->update();
//    d_fissionsource->setup_outer(1.0 / keff);
//    State::moments_type fd_old(d_fissionsource->density());
//
//    // Update density.
//    d_fissionsource->update();

    // Compute keff.
    const State::moments_type &fd = d_fissionsource->density();
    keff *= norm(fd, "L1") / norm(fd_old, "L1");

    // Compute error in fission density.
    error = norm_residual(fd, fd_old, "L1");
    if (d_print_level > 1 && iteration % d_print_interval == 0)
    {
      printf("PI Iter: %3i  Error: %12.9f  keff: %12.9f  cmfd: %12.9f \n",
             iteration, error, keff, keff_c);
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

  d_state->set_eigenvalue(keff);
  std::cout << "CMFD done." << std::endl;
}

//----------------------------------------------------------------------------//
template <class D>
double EigenCMFD<D>::cmfd_update()
{
  // Get the solver explicitly as MG-CMFD
  typedef MGSolverCMFD<D> Solver_T;
  Solver_T *solver = dynamic_cast<Solver_T*>(&(*d_mg_solver->solver()));
  SP_mesh coarse_mesh = solver->coarse_mesh();

  // Homogenize the material
  Homogenize H(d_material, Homogenize::PHI_D);
  SP_material cmat = H.homogenize(d_state, d_mesh, "COARSEMESH");
  const vec2_dbl &phi = H.coarse_mesh_flux();

  // Create loss matrix
  typedef typename CMFDLossOperator<D>::SP_lossoperator SP_L;
  SP_L L(new CMFDLossOperator<D>(d_input,
                                 cmat,
                                 solver->coarse_mesh(),
                                 solver->tally(),
                                 false, // not multiplying
                                 d_adjoint));
  L->construct(phi);

  // Create gain matrix
  SP_matrix F(new DiffusionGainOperator(d_input,
                                        cmat,
                                        coarse_mesh,
                                        d_adjoint));

  // Construct source
  callow::Vector x(d_number_groups * coarse_mesh->number_cells(), 0.0);
  const vec_int &cmap = d_mesh->mesh_map("COARSEMESH");
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t ci = cmap[i] + g * coarse_mesh->number_cells();
      double v_ratio = d_mesh->volume(i) / coarse_mesh->volume(cmap[i]);
      x[ci] += d_state->phi(g)[i] * v_ratio;
    }
  }
  x.scale(1.0/x.norm());
  callow::Vector x0(x);

  // Solve the eigenvalue problem
  d_eigensolver->set_operators(F, L, d_eigen_db);
  d_eigensolver->solve(x, x0);

  // Scale the fluxes
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t ci = cmap[i] + g * coarse_mesh->number_cells();
      d_state->phi(g)[i] *= (d_omega * x[ci] / x0[ci] + (1-d_omega));
    }
  }

  return d_eigensolver->eigenvalue();
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class EigenCMFD<_1D>;
template class EigenCMFD<_2D>;
template class EigenCMFD<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of EigenCMFD.cc
//----------------------------------------------------------------------------//
