//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Richardson.cc
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  Richardson class definition.
 */
//---------------------------------------------------------------------------//

#include "Richardson.hh"

namespace callow
{

//---------------------------------------------------------------------------//
Richardson::Richardson(const double  atol,
                       const double  rtol,
                       const int     maxit,
                       const double  omega)
  : LinearSolver(atol, rtol, maxit, "richardson")
  , d_omega(omega)
{
  /* ... */
}

} // end namespace callow

