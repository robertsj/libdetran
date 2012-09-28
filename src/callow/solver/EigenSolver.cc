//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EigenSolver.cc
 *  @brief  EigenSolver member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 26, 2012
 */
//---------------------------------------------------------------------------//

#include "EigenSolver.hh"

namespace callow
{

EigenSolver::EigenSolver(const double    tol,
                         const int       maxit,
                         std::string     name)
  : d_tolerance(tol)
  , d_maximum_iterations(maxit)
  , d_name(name)
  , d_residual_norm(maxit + 1, 0)
  , d_number_iterations(0)
  , d_monitor_level(2)
{
  Require(d_tolerance >= 0.0);
  Require(d_maximum_iterations > 0);
}

} // end namespace callow

//---------------------------------------------------------------------------//
//              end of file EigenSolver.cc
//---------------------------------------------------------------------------//
