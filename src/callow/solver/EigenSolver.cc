//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenSolver.cc
 *  @brief EigenSolver member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenSolver.hh"
#include "LinearSolverCreator.hh"
#include "preconditioner/PCILU0.hh"
#define W() std::cout << "line " << __LINE__ << " on file " << __FILE__ << std::endl;

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
  d_A->number_rows();
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
int EigenSolver::solve(Vector &x, Vector &x0)
{
  Insist(d_A, "The operator Ad cannot be null");
  d_A->number_rows();
  Require(x.size() == d_A->number_rows());
  if (x0.size())
    Require(x0.size() == d_A->number_rows());
  d_status = MAXIT;
  solve_impl(x, x0);
  if (d_status ==  MAXIT && d_monitor_level > 0)
  {
    printf("*** %s did not converge within the maximum number of iterations\n",
           d_name.c_str());
  }
  return d_status;
}

//----------------------------------------------------------------------------//
bool EigenSolver::monitor(int it, double l, double r)
{
  // record the iteration and residual norm
  d_number_iterations = it;
  d_residual_norm[it] = r;
  // echo the residual
  if (d_monitor_level > 1)
  {
    printf("iteration: %5i  eigenvalue: %12.8e    residual: %12.8e \n",
           it, l, r);
  }
  // send a signal
  if (it == d_maximum_iterations)
  {
    d_status = MAXIT;
    return true;
  }
  if (r < d_tolerance)
  {
    if (d_monitor_level > 0)
    {
      printf("*** %s converged in %5i iterations with a residual of %12.8e \n",
             d_name.c_str(), it, r);
    }
    d_status = SUCCESS;
    return true;
  }
  return false;
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
