//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionEigensolver.cc
 * \author robertsj
 * \date   Jul 26, 2012
 * \brief  DiffusionEigensolver class definition.
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_SLEPC

// Diffusion
#include "DiffusionEigensolver.hh"

// System
#include <string>
#include <cmath>
#include <iostream>

namespace detran_diffusion
{

DiffusionEigensolver::DiffusionEigensolver(SP_input    input,
                                           SP_material material,
                                           SP_mesh     mesh,
                                           SP_state    state)
  : d_input(input)
  , d_state(state)
  , d_tolerance(1e-5)
  , d_max_iters(10)
{
  // Preconditions
  Require(d_input);
  Require(material);
  Require(mesh);

  // Create loss operator
  d_M = new LossOperator(input, material, mesh);

  // Create gains operator
  d_F = new GainOperator(input, material, mesh);

  // Number of groups
  d_number_groups =  material->number_groups();

  // Group size
  d_group_size = mesh->number_cells();

  // System size
  d_size = mesh->number_cells() * material->number_groups();

  PetscErrorCode ierr;

  // Create unknown vector
  ierr = VecCreateSeq(PETSC_COMM_SELF, d_size, &d_x);
  Insist(!ierr, "Error creating Vec.");

  // Create the eigensolver
  ierr = EPSCreate(PETSC_COMM_SELF, &d_solver);
  Insist(!ierr, "Error creating EPS context.");

  // Set the operator.
  ierr = EPSSetOperators(d_solver, d_F->get_operator(), d_M->get_operator());
  Insist(!ierr, "Error setting EPS operator.");

  // Set the problem type.
  ierr = EPSSetProblemType(d_solver, EPS_GNHEP);
  Insist(!ierr, "Error setting EPS problem type.");

  // Set the solver type to krylovschur and one largest eigenvalue.
  ierr = EPSSetType(d_solver, EPSKRYLOVSCHUR);
  Insist(!ierr, "Error defaulting EPS to EPSKRYLOVSCHUR.");
  ierr = EPSSetWhichEigenpairs(d_solver, EPS_LARGEST_MAGNITUDE);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_solver);
  Insist(!ierr, "Error setting EPS from options.");

  // Set convergence criteria
  if (input->check("diffusion_max_iters"))
    d_max_iters = input->get<int>("diffusion_max_iters");
  if (input->check("diffusion_tolerance"))
    d_tolerance = input->get<double>("diffusion_tolerance");

  // Set the convergence criteria
  ierr = EPSSetTolerances(d_solver, d_tolerance, d_max_iters);
  Insist(!ierr, "Error setting EPS tolerances.");

  // Postconditions
  Ensure(d_M);
  Ensure(d_F);
}

void DiffusionEigensolver::solve()
{
  using std::cout;
  using std::endl;

  cout << "Starting DIFFUSION EIGENSOLVE." << endl;

  // Temporaries
  double keff_imag;
  Vec    x_imag;
  VecDuplicate(d_x, &x_imag);
  int    ierr;

  // Solve the eigenproblem
  ierr = EPSSolve(d_solver);
  Insist(!ierr, "Error solving current eigenvalue problem.");

  // Get the number of iterations.
  int numit = 0;
  ierr = EPSGetIterationNumber(d_solver, &numit);
  Insist(!ierr, "Error getting iteration count.");
  std::cout << "EIGEN: Number of iterations =" << numit << std::endl;

  // Check the number of converged.
  int numconv = 0;
  ierr = EPSGetConverged(d_solver, &numconv);
  Insist(!ierr, "Error getting iteration count.");
  if (numconv == 0)
  {
    std::cout << "No converged eigenvalues!" << std::endl;
  }
  else
  {

    // Get the dominant mode.
    ierr = EPSGetEigenpair(d_solver, 0, &d_keff, &keff_imag, d_x, x_imag);
    Insist(!ierr, "Error getting eigenpair.");

    // Scale the result by its sum.  This points it in the positive
    // direction and gives it an L1 normalization.
    double sum;
    VecSum(d_x, &sum);
    VecScale(d_x, 1.0/sum);

    // Fill the state.
    fill_state();
  }

  // Free temporary
  VecDestroy(&x_imag);

  //
  std::cout << "EIGEN done. keff = " << d_keff <<  std::endl;

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

void DiffusionEigensolver::fill_state()
{
  d_state->set_eigenvalue(d_keff);

  // Get the Vec array
  double *x_a;
  VecGetArray(d_x, &x_a);

  // Fill the state
  for (int g = 0; g < d_number_groups; g++)
  {
    for (int i = 0; i < d_group_size; i++, x_a++)
    {
      d_state->phi(g)[i] = *x_a;
    }
  }

  VecRestoreArray(d_x, &x_a);

}


} // end namespace detran_diffusion

#endif

//---------------------------------------------------------------------------//
//              end of file DiffusionEigensolver.cc
//---------------------------------------------------------------------------//
