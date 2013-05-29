//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  JacobiDavidson.cc
 *  @brief JacobiDavidson member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "JacobiDavidson.hh"

namespace callow
{

//----------------------------------------------------------------------------//
JacobiDavidson::JacobiDavidson(const double    tol,
                               const int       maxit,
                               const int       subspace_size)
  : Base(tol, maxit, "power")
  , d_subspace_size(subspace_size)
  , d_theta(0.0)
{
  /* ... */
}

//----------------------------------------------------------------------------//
void JacobiDavidson::set_operators(SP_matrix A,
                                   SP_matrix B,
                                   SP_db     db)
{
  EigenSolver::set_operators(A, B, db);
  // Create the default preconditioner
  d_P = new JacobiDavidsonDefaultP(A, B, d_solver, this);
}

//----------------------------------------------------------------------------//
void JacobiDavidson::set_preconditioner(SP_matrix P)
{
  Require(P);
  d_P = P;
}

//----------------------------------------------------------------------------//
void JacobiDavidson::solve_impl(Vector &x, Vector &x0)
{

}

} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file JacobiDavidson.cc
//----------------------------------------------------------------------------//



