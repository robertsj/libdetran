//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  LinearSolver.i.hh
 *  @brief LinearSolver inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_LINEARSOLVER_I_HH_
#define callow_LINEARSOLVER_I_HH_

namespace callow
{

//----------------------------------------------------------------------------//
inline int LinearSolver::solve(const Vector &b, Vector &x)
{
  Require(x.size() == b.size());
  Require(x.size() == d_A->number_rows());
  // Resize the norm
  d_residual.resize(d_maximum_iterations+1, 0.0);
  d_status = RUNNING;
  solve_impl(b, x);
  if (d_status ==  MAXIT && d_monitor_level > 0)
  {
     printf("*** %s did not converge within the maximum number of iterations\n",
            d_name.c_str());
  }
  // Resize the norm
  d_residual.resize(d_number_iterations+1);
  return d_status;
}

} // end namespace callow

#endif // callow_LINEARSOLVER_I_HH_

//----------------------------------------------------------------------------//
//              end of file LinearSolver.i.hh
//----------------------------------------------------------------------------//
