//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobi.cc
 *  @brief Jacobi class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Jacobi.hh"

namespace callow
{

//----------------------------------------------------------------------------//
Jacobi::Jacobi(const double  atol,
               const double  rtol,
               const int     maxit,
               const double  omega,
               bool          successive_norm)
  : LinearSolver(atol, rtol, maxit, "jacobi")
  , d_successive_norm(successive_norm)
  , d_omega(omega)
{
  /* ... */
}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file Jacobi.cc
//----------------------------------------------------------------------------//
