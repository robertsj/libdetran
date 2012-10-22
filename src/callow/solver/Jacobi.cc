//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Jacobi.cc
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Jacobi class definition.
 */
//---------------------------------------------------------------------------//

#include "Jacobi.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

Jacobi::Jacobi(const double  atol,
               const double  rtol,
               const int     maxit,
               bool successive_norm)
  : LinearSolver(atol, rtol, maxit, "jacobi")
  , d_successive_norm(successive_norm)
{
  /* ... */
}

} // end namespace callow

