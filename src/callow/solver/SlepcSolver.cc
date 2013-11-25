//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SlepcSolver.cc
 *  @brief SlepcSolver member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "callow/callow_config.hh"
#include "SlepcSolver.hh"

namespace callow
{

#ifdef CALLOW_ENABLE_SLEPC

//----------------------------------------------------------------------------//
SlepcSolver::SlepcSolver(const std::string &eps_type,
                         const double       tol,
                         const int          maxit,
                         const int          number_values)
  : EigenSolver(tol, maxit, "solver_slepc")
  , d_eps_type(eps_type)
{
  PetscErrorCode ierr;

  // Create the context.
  ierr = EPSCreate(PETSC_COMM_SELF, &d_slepc_solver);
  Insist(!ierr, "Error creating EPS context.");

  // Set the solver type and get just the largest eigenvalue.
  ierr = EPSSetWhichEigenpairs(d_slepc_solver, EPS_LARGEST_MAGNITUDE);
  Insist(!ierr, "Error selecting EPS eigenpairs.");

  // Set the tolerances
  ierr = EPSSetTolerances(d_slepc_solver, d_tolerance, d_maximum_iterations);
  Insist(!ierr, "Error setting EPS tolerances.");
}

//----------------------------------------------------------------------------//
SlepcSolver::~SlepcSolver()
{
  EPSDestroy(&d_slepc_solver);
}

//----------------------------------------------------------------------------//
void SlepcSolver::solve_impl(Vector &x, Vector &x0)
{
  double lambda;
  double lambda_imag;

  // Solve the eigenproblem
  PetscErrorCode ierr = EPSSolve(d_slepc_solver);
  Insist(!ierr, "Error solving eigenvalue problem.");

  // Extract the number of iterations
  ierr = EPSGetIterationNumber(d_slepc_solver, &d_number_iterations);
  Insist(!ierr, "Error getting iteration count.");

  // Get the dominant mode.
  ierr = EPSGetEigenpair(d_slepc_solver, 0, &lambda, &lambda_imag,
                         x.petsc_vector(), x0.petsc_vector());
  Insist(!ierr, "Error getting eigenpair.");

  // Scale the result by its sum.  This points it in the positive
  // direction and gives the L1 normalization we're using in PI.
  PetscScalar sign = 1;
  if (x[0] < 0) sign = -1;
  x.scale(sign / x.norm(L1));

  // Store the eigenvalue
  d_lambda = lambda;
}

//----------------------------------------------------------------------------//
void SlepcSolver::set_operators(SP_matrix    A,
                                SP_matrix    B,
                                SP_db        db)
{
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
    ierr = EPSSetOperators(d_slepc_solver,
                           d_A->petsc_matrix(), d_B->petsc_matrix());
    Insist(!ierr, "Error setting EPS operator.");
  }

  // Set the solver type
  set_eps_type(d_eps_type);

  // Then allow for user choice.
  ierr = EPSSetFromOptions(d_slepc_solver);
  Insist(!ierr, "Error setting EPS from options.");
}

//----------------------------------------------------------------------------//
void SlepcSolver::set_preconditioner(SP_preconditioner  P,
                                     const int          side)
{
  Require(P);
  d_P = P;
  ST st;
  PetscErrorCode ierr;
  ierr = EPSGetST(d_slepc_solver, &st);
  //ierr = STSetType(st, STPRECOND);
  d_P->set_slepc_st(st);
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void SlepcSolver::set_preconditioner_matrix(SP_matrix  P,
                                            const int  side)
{
  Require(P);
  ST st;
  PetscErrorCode ierr = EPSGetST(d_slepc_solver, &st);
  // set the user-given matrix as the preconditioning matrix
  ierr = STPrecondSetMatForPC(st, P->petsc_matrix());
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void SlepcSolver::set_eps_type(const std::string &eps_type)
{
  d_eps_type = eps_type;
  // Set the solver type
  PetscErrorCode ierr = EPSSetType(d_slepc_solver, d_eps_type.c_str());
  Insist(!ierr, "Error creating EPS of type " + d_eps_type);
}

#else

SlepcSolver::SlepcSolver(const std::string&,const double,const int,const int)
  : d_eps_type(""), d_slepc_solver(NULL) {}
SlepcSolver::~SlepcSolver(){}
void SlepcSolver::solve_impl(Vector &, Vector &){}
void SlepcSolver::set_operators(SP_matrix, SP_matrix, SP_db){}
void SlepcSolver::set_preconditioner(SP_preconditioner,const int){}
void SlepcSolver::set_preconditioner_matrix(SP_matrix, const int){}
void SlepcSolver::set_eps_type(const std::string &){}

#endif // CALLOW_ENABLE_SLEPC

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file SlepcSolver.cc
//----------------------------------------------------------------------------//
