//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LinearSolver.cc
 * \brief  LinearSolver member definitions
 * \author Jeremy Roberts
 * \date   Sep 26, 2012
 */
//---------------------------------------------------------------------------//

#include "LinearSolver.hh"

namespace callow
{

//-------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//-------------------------------------------------------------------------//

LinearSolver::LinearSolver(const double atol,
                           const double rtol,
                           const int    maxit,
                           std::string  name)
  : d_absolute_tolerance(atol)
  , d_relative_tolerance(rtol)
  , d_maximum_iterations(maxit)
  , d_residual(maxit + 1, 0)
  , d_number_iterations(0)
  , d_monitor_level(2)
  , d_monitor_diverge(true)
  , d_norm_type(L2)
  , d_name(name)
  , d_pc_side(NONE)
{
  Require(d_absolute_tolerance >= 0.0);
  Require(d_relative_tolerance >= 0.0);
  Require(d_maximum_iterations >  0);
}

//-------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//-------------------------------------------------------------------------//

void LinearSolver::
set_operators(SP_matrix A,
              SP_preconditioner P,
              const int side)
{
  Require(A);
  d_A = A;
  Ensure(d_A->number_rows() == d_A->number_columns());
  if (P) d_P = P;
  d_pc_side = side;
}

void LinearSolver::
set_tolerances(const double atol, const double rtol, const int maxit)
{
  d_absolute_tolerance = atol;
  d_relative_tolerance = rtol;
  d_maximum_iterations = maxit;
  Require(d_absolute_tolerance > 0.0);
  Require(d_relative_tolerance > 0.0);
  Require(d_maximum_iterations >= 0);
}

// print out iteration and residual for initial
bool LinearSolver::monitor_init(double r)
{
  d_residual[0] = r;
  if (d_monitor_level > 1)
    printf("iteration: %5i    residual: %12.8e \n", 0, r);
  if (r < d_absolute_tolerance)
  {
    if (d_monitor_level > 0)
    {
      printf("*** %s converged in %5i iterations with a residual of %12.8e \n",
             d_name.c_str(), 0, r );
    }
    d_status = SUCCESS;
    return true;
  }
  return false;
}

// print out iteration and residual
bool LinearSolver::monitor(int it, double r)
{
  d_number_iterations = it;
  d_residual[it] = r;
  if (d_monitor_level > 1)
    printf("iteration: %5i    residual: %12.8e \n", it, r);
 // Assert(it > 0);
  if (r < std::max(d_relative_tolerance * d_residual[0],
                   d_absolute_tolerance))
  {
    if (d_monitor_level)
    {
      printf("*** %s converged in %5i iterations with a residual of %12.8e \n",
             d_name.c_str(), it, r );
    }
    d_status = SUCCESS;
    return true;
  }
  else if (d_monitor_diverge and it >  1 and r - d_residual[it - 1] > 0.0)
  {
    if (d_monitor_level) printf("*** %s diverged \n", d_name.c_str());
    d_status = DIVERGE;
    return true;
  }
  return false;
}


} // end namespace callow

//---------------------------------------------------------------------------//
//              end of file LinearSolver.cc
//---------------------------------------------------------------------------//
