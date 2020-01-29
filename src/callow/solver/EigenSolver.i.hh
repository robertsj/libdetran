//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenSolver.i.hh
 *  @brief EigenSolver inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_EIGENSOLVER_I_HH_
#define callow_EIGENSOLVER_I_HH_

#include "LinearSolverCreator.hh"

namespace callow
{

//----------------------------------------------------------------------------//
inline int EigenSolver::solve(Vector &x, Vector &x0)
{
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
inline bool EigenSolver::monitor(int it, double l, double r)
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

} // end namespace callow

#endif /* callow_EIGENSOLVER_I_HH_ */


//----------------------------------------------------------------------------//
//              end of file EigenSolver.i.hh
//----------------------------------------------------------------------------//
