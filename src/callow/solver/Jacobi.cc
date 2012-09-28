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
               const int     maxit)
  : LinearSolver(atol, rtol, maxit, "jacobi")
{
  /* ... */
}

} // end namespace callow

