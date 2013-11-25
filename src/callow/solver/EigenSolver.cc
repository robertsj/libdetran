//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenSolver.cc
 *  @brief EigenSolver member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenSolver.hh"
#include "preconditioner/PCILU0.hh"

namespace callow
{

//----------------------------------------------------------------------------//
EigenSolver::EigenSolver(const double    tol,
                         const int       maxit,
                         std::string     name)
  : d_tolerance(tol)
  , d_maximum_iterations(maxit)
  , d_name(name)
  , d_residual_norm(maxit + 1, 0)
  , d_number_iterations(0)
  , d_monitor_level(2)
  , d_lambda(0.0)
  , d_status(RUNNING)
{
  Require(d_tolerance >= 0.0);
  Require(d_maximum_iterations > 0);
}

//----------------------------------------------------------------------------//
void EigenSolver::set_operators(SP_matrix    A,
                                SP_matrix    B,
                                SP_db        db)
{
  Insist(A, "The operator A cannot be null");
  d_A = A;
  Ensure(d_A->number_rows() == d_A->number_columns());

  // Setup linear system if this is a generalized eigenproblem
  if (B)
  {
    d_B = B;
    Ensure(d_B->number_rows() == d_B->number_columns());
    Ensure(d_B->number_rows() == d_A->number_columns());
    // Create linear solver.  Defaults can be changed by
    // extracting the solver and resetting.  The same goes
    // for the preconditioner, which is null by default.
    d_solver = LinearSolverCreator::Create(db);
    d_solver->set_operators(d_B, db);
  }
}

//----------------------------------------------------------------------------//
void EigenSolver::set_preconditioner(SP_preconditioner P, const int side)
{
  Require(P);
  // By default, the preconditioner is passed on to the linear solver
  d_solver->set_preconditioner(P, side);
}

//----------------------------------------------------------------------------//
void EigenSolver::set_preconditioner_matrix(SP_matrix P, const int side)
{
  Require(P);
  SP_preconditioner PC(new PCILU0(P));
  d_solver->set_preconditioner(PC, side);
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file EigenSolver.cc
//----------------------------------------------------------------------------//
