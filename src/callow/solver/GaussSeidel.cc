//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussSeidel.cc
 * \author robertsj
 * \date   Sep 14, 2012
 * \brief  GaussSeidel class definition.
 */
//---------------------------------------------------------------------------//

#include "GaussSeidel.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

GaussSeidel::GaussSeidel(const double  atol,
                         const double  rtol,
                         const int     maxit,
                         const double  omega,
                         bool          successive_norm)
  : LinearSolver(atol, rtol, maxit, "gauss-seidel")
  , d_omega(omega)
  , d_successive_norm(successive_norm)
{
  /* ... */
}

} // end namespace callow



