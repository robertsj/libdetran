//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SlepcSolver.cc
 *  @author robertsj
 *  @date   Sep 25, 2012
 *  @brief  SlepcSolver class definition.
 */
//---------------------------------------------------------------------------//

#include "callow/callow_config.hh"

#ifdef CALLOW_ENABLE_SLEPC

#include "SlepcSolver.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

SlepcSolver::SlepcSolver(const double  tol,
                         const int     maxit)
  : EigenSolver<PetscScalar>(tol, maxit, "solver_slepc")
{
  PetscErrorCode ierr;
  std::cout << " start slepcsolver... " << std::endl;
  // Create the context.
  ierr = EPSCreate(PETSC_COMM_SELF, &d_slepc_solver);
  Insist(!ierr, "Error creating EPS context.");

  // Set the solver type to krylovschur and one largest eigenvalue.
  ierr = EPSSetType(d_slepc_solver, EPSKRYLOVSCHUR);
  Insist(!ierr, "Error defaulting EPS to EPSKRYLOVSCHUR.");
  ierr = EPSSetWhichEigenpairs(d_slepc_solver, EPS_LARGEST_MAGNITUDE);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Set the tolerances
  ierr = EPSSetTolerances(d_slepc_solver, d_tolerance, d_maximum_iterations);
  Insist(!ierr, "Error setting EPS tolerances.");

}

SlepcSolver::~SlepcSolver()
{
  std::cout << "*** DESTROYING SLEPC ***" << std::endl;
  EPSDestroy(&d_slepc_solver);
}

void SlepcSolver::set_operators(SP_matrix   A,
                                SP_matrix   B,
                                std::string type)
{
  // Preconditions
  Insist(A, "First operator cannot be null");

  // Store the operators
  d_A = A;
  if (B) d_B = B;

  // Set the operators.
  PetscErrorCode ierr;
  if (!B)
  {
    // standard Ax=Lx, nonhermitian evp
    ierr = EPSSetProblemType(d_slepc_solver, EPS_NHEP);
    Insist(!ierr, "Error setting EPS problem type.");
    ierr = EPSSetOperators(d_slepc_solver, d_A->petsc_matrix(), PETSC_NULL);
    Insist(!ierr, "Error setting EPS operator.");
  }
  else
  {
    // generalized Ax=LBx, nonhermitian evp
    ierr = EPSSetProblemType(d_slepc_solver, EPS_GNHEP);
    Insist(!ierr, "Error setting EPS problem type.");
    ierr = EPSSetOperators(d_slepc_solver, d_A->petsc_matrix(), d_B->petsc_matrix());
    Insist(!ierr, "Error setting EPS operator.");
  }

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_slepc_solver);
  Insist(!ierr, "Error setting EPS from options.");

}

} // end namespace callow

#endif // CALLOW_ENABLE_SLEPC
