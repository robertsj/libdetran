//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   DiffusionFixedSourceSolve.cc
 * \brief  DiffusionFixedSourceSolve member definitions
 * \author Jeremy Roberts
 * \date   Sep 11, 2012
 */
//---------------------------------------------------------------------------//

#include "DiffusionFixedSourceSolver.hh"

namespace detran
{

template <class D>
DiffusionFixedSourceSolver<D>::DiffusionFixedSourceSolver(SP_input input,
                                                          SP_material material,
                                                          SP_mesh mesh,
                                                          SP_state state)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_state(state)
  , d_boundary(new BoundaryDiffusion<D>(input, mesh))
  , d_fixed_type(FIXED)
  , d_maximum_iterations(100)
  , d_maximum_fission_iterations(20)
  , d_tolerance(1e-6)
  , d_fission_tolerance(1e-6)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_state);
  Require(d_boundary);

  d_problem_size = d_mesh->number_cells() * d_material->number_groups();

  // Select solver type
  if (d_input->check("diffusion_fixed_type"))
  {
    d_fixed_type = d_input->template get<int>("diffusion_fixed_type");
    Ensure(d_fixed_type < END_FIXED_TYPES);
  }

  // Create operators and vectors
  if (d_fixed_type == MULTIPLY)
    d_M = new DiffusionLossOperator(d_input, d_material, d_mesh, true, 1.0);
  else
    d_M =  new DiffusionLossOperator(d_input, d_material, d_mesh, false);
  if (d_fixed_type == ITERATE)
    d_F = new DiffusionGainOperator(d_input, d_material, d_mesh);
  d_x = new Vector(d_problem_size, 0.0);
  d_b = new Vector(d_problem_size, 0.0);

  // Create solver
  PetscErrorCode ierr;
  ierr = KSPCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating KSP object.");

  // Set the operators.  We use the same operator for the PC.
  KSPSetOperators(d_solver, d_M->A(), d_M->A(), SAME_NONZERO_PATTERN);

  // Set tolerances.
  ierr = KSPSetTolerances(d_solver,
                          d_tolerance,   // relative tolerance
                          PETSC_DEFAULT, // absolute tolerance
                          PETSC_DEFAULT, // divergence tolerance
                          d_maximum_iterations);

  // Allow for command line flags.
  ierr = KSPSetFromOptions(d_solver);

  // Postconditions
  Ensure(!ierr);
}

template <class D>
void DiffusionFixedSourceSolver<D>::solve()
{
  // Call the appropriate solver
  if (d_fixed_type == FIXED or d_fixed_type == MULTIPLY)
    solve_fixed();
  else if (d_fixed_type == ITERATE)
    solve_iterate();
  else
    THROW("Unsupported diffusion fixed source type");

  // Fill the state with the converged flux

  // Fill the boundary

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

// Direct solve
template <class D>
void DiffusionFixedSourceSolver<D>::solve_fixed()
{
  PetscErrorCode ierr = KSPSolve(d_solver, d_b, d_x);
  Insist(!ierr, "Error in KSPSolve.");
}

// Fission source iterations
template <class D>
void DiffusionFixedSourceSolver<D>::solve_iterate()
{
  PetscErrorCode ierr = KSPSolve(d_solver, d_b, d_x);
  Insist(!ierr, "Error in KSPSolve.");
}



} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file DiffusionFixedSourceSolve.cc
//---------------------------------------------------------------------------//
