//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenCMFD.cc
 *  @brief EigenCMFD class definition
 *  @note  Copyright(C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "solvers/eigen/EigenCMFD.hh"
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
}

//----------------------------------------------------------------------------//
template <class D>
void EigenCMFD<D>::solve()
{
  using detran_utilities::norm_residual;
  using detran_utilities::norm;
  using detran_utilities::vec_scale;

  typedef MGSolverCMFD<D> Solver_T;
  Solver_T *solver = dynamic_cast<Solver_T*>(&(*d_mg_solver->solver()));
  Insist(solver, "EigenCMFD needs CMFD as the MG problem to get operator.");

  std::cout << "Starting CMFD." << std::endl;

  double keff = 1.0, keff_old;

  // Initialize density with normalize (L1) guess
  d_fissionsource->initialize();

  double error_k, error_fd;
  int iteration;

  for (iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {
    // Do one energy sweep.
    d_fissionsource->setup_outer(1.0);
    solver->energy_sweep();

    // CMFD update
    keff_old = keff;
    keff = cmfd_update();

    // Compute errors
    State::moments_type fd_old(d_fissionsource->density());
    d_fissionsource->update();
    State::moments_type &fd = d_fissionsource->density();
    vec_scale(fd, 1.0 / norm(fd, "L1"));
    error_fd = norm_residual(fd, fd_old, "Linf");
    error_k  = (keff - keff_old) / keff_old;

    // Monitor
    if (d_print_level > 1 && iteration % d_print_interval == 0)
    {
      printf("CMFD ITER:  %3i  keff: %12.9f  err_k: %12.5e  err_fd: %12.5e \n",
             iteration, keff, error_k, error_fd);
    }

    if (error_fd < d_tolerance) break;

  } // eigensolver loop

  if (d_print_level > 0)
  {
    std::string line = "";
    for (int i = 0; i < 80; ++i) line += "*";
    line += "\n";
    printf(line.c_str());
    printf("CMFD FINAL: %3i  keff: %12.9f  err_k: %12.5e  err_fd: %12.5e \n",
           iteration, keff, error_k, error_fd);
    printf(line.c_str());
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

  //compute_current();
  Homogenize H(d_material, Homogenize::PHI_D);


  SP_material cmat = H.homogenize(d_state, d_mesh, "COARSEMESH");
  const vec2_dbl &phi = H.coarse_mesh_flux();

  // Create loss matrix
  if (!d_loss)
  {
    d_loss = new CMFDLossOperator<D>(d_input,
                                     cmat,
                                     solver->coarse_mesh(),
                                     solver->tally(),
                                     false, // not multiplying
                                     d_adjoint);
    d_loss->construct(phi, 1.0, cmat, true);
  }
  else
  {
    d_loss->construct(phi, 1.0, cmat, false);
  }

  // Create gain matrix
  if (!d_gain)
  {
    d_gain = new DiffusionGainOperator(d_input,
                                       cmat,
                                       coarse_mesh,
                                       d_adjoint);
  }
  else
  {
    d_gain->construct(cmat);
  }

  // Construct coarse mesh fluxes
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
  x.scale(1.0 / x.norm());
  callow::Vector x0(x);
  callow::Vector xI(x); // callow/slepc overwrites initial guess

  // Solve the eigenvalue problem
  d_eigensolver = callow::EigenSolverCreator::Create(d_eigen_db);
  d_eigensolver->set_operators(d_gain, d_loss, d_eigen_db);
  d_eigensolver->solve(x, xI);
  x.scale(1.0 / x.norm());
  static int count = 0;
  //x.print_matlab("X"+AsString(count)+".out");

  ++count;

  // Scale the fluxes
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      size_t ci = cmap[i] + g * coarse_mesh->number_cells();
      d_state->phi(g)[i] *= (d_omega * x[ci] / x0[ci] + (1.0-d_omega));
    }
  }

  return d_eigensolver->eigenvalue();
}

//----------------------------------------------------------------------------//
template <class D>
void EigenCMFD<D>::compute_current()
{

  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      double J[3] = {0.0, 0.0, 0.0};
      for (size_t o = 0; o < d_mg_solver->quadrature()->number_octants(); ++o)
      {
        for (size_t a = 0; a < d_mg_solver->quadrature()->number_angles_octant(); ++a)
        {
          double w = d_mg_solver->quadrature()->weight(a);
          double psi = d_state->psi(g, o, a)[i];
          for (size_t d = 0; d < D::dimension; ++d)
          {
            J[d] += w * psi * d_mg_solver->quadrature()->cosines(d)[a];
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

template class EigenCMFD<_1D>;
template class EigenCMFD<_2D>;
template class EigenCMFD<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of EigenCMFD.cc
//----------------------------------------------------------------------------//
