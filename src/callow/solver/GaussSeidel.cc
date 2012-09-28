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
                         const int     maxit)
  : LinearSolver(atol, rtol, maxit, "gauss-seidel")
{
  /* ... */
}

} // end namespace callow



