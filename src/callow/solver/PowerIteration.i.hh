//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  PowerIteration.i.hh
 *  @brief PowerIteration inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef callow_POWERITERATION_I_HH_
#define callow_POWERITERATION_I_HH_

#include <cmath>
#include <iostream>

namespace callow
{

inline void PowerIteration::solve_impl(Vector &x, Vector &x0)
{
  // Initialize guess if not present
  if (!x0.size()) x0.resize(x.size(), 1.0);

  // Eigenvalue and error.
  double lambda_0   = 0.0;
  double lambda     = 0.0;
  double lambda_err = 0.0;

  // Temporary vector
  Vector temp(x0);

  // Eigenvector errors.  We'll estimate the dominance ratio using
  // a few iterations, and then estimate how many iterations until
  // the error has died out such that the tolerance is met.
  double x_err_0  = 1.0;
  double x_err    = 1.0;

  // Normalize guess to unity
  lambda = x0.norm(L1);
  x0.scale(1.0 / lambda);

  // Do power iterations
  for (int i = 0; i < d_maximum_iterations; ++i)
  {

    // Save previous value.
    lambda_0 = lambda;
    x_err_0  = x_err;

    // compute x <-- A*x0
    d_A->multiply(x0, x);

    // If B exists, solve B*x = A*x0
    if (d_B)
    {
      d_solver->solve(x, temp);
      x.copy(temp);
    }

    // Compute the eigenvalue and update eigenvalue error.
    lambda     = x.norm(L1);
    lambda_err = std::fabs(lambda - lambda_0);

    // Compute residual
    x0.scale(-lambda_0);
    x0.add(x);
    x_err = x0.norm(L1);

    // Normalize the vector and save.
    x.scale(1.0 / lambda);
    x0.copy(x);

    // Check for convergence.
    if (monitor(i, lambda, x_err)) break;

  }

  /// Store the eigenvalue
  d_lambda = lambda;

}

} // end namespace callow

#endif // callow_POWERITERATION_I_HH_

//---------------------------------------------------------------------------//
//              end of file PowerIteration.i.hh
//---------------------------------------------------------------------------//
