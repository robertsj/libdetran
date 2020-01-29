//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobi.i.hh
 *  @brief Jacobi inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_JACOBI_I_HH_
#define callow_JACOBI_I_HH_

#include "callow/matrix/Matrix.hh"
#include <cmath>
#include <cstdio>
#ifdef DETRAN_ENABLE_OPENMP
#include <omp.h>
#endif

namespace callow
{

//----------------------------------------------------------------------------//
// SOLVE
//----------------------------------------------------------------------------//

inline void Jacobi::solve_impl(const Vector &b, Vector &x)
{
  Insist(dynamic_cast<Matrix*>(d_A.bp()),
    "Need an explicit matrix for use with Jacobi iteration");
  Matrix::SP_matrix A = d_A;

  // temporary storage and pointers for swapping
  Vector temp(x.size(), 0.0);
  Vector* x0 = &x;
  Vector* x1 = &temp;
  Vector* swap;

  // iteration count
  int &iteration = d_number_iterations;

  // compute initial residual Ax - b and its norm
  A->multiply((*x0), (*x1));
  double r = x1->norm_residual(b, L2);
  if (monitor_init(r))
  {
    //return;
  }

  // perform iterations
  for (iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {

    //----------------------------------------------------//
    // compute X1 <-- -inv(D)*(L+U)*X0 + inv(D)*b
    //----------------------------------------------------//

    double* a = A->values();
    #pragma omp for
    for (int i = 0; i < A->number_rows(); ++i)
    {
      double v = 0;
      int p = A->start(i);
      int d = A->diagonal(i);
      // L * X0
      for (; p < d; ++p)
        v += a[p] * (*x0)[A->column(p)];
      ++p; // skip diagonal
      // U * X0
      for (; p < A->end(i); ++p)
        v += a[p] * (*x0)[A->column(p)];
      // weighted result
      (*x1)[i] = d_omega * (b[i] - v) / a[d] + (1.0 - d_omega) * (*x0)[i];
    }
    a = 0; // nullify pointer

    //----------------------------------------------------//
    // compute residual norm
    //----------------------------------------------------//

    if (!d_successive_norm)
    {
      // Compute the residual and put it into x0
      A->multiply(*x1, *x0);
      r = x0->norm_residual(b, d_norm_type);
    }
    else
    {
      // compare x1 with x0
      r = x1->norm_residual(*x0, d_norm_type);
    }

    //----------------------------------------------------//
    // swap pointers
    //----------------------------------------------------//
    swap = x0;
    x0   = x1;
    x1   = swap;

    //----------------------------------------------------//
    // check convergence
    //----------------------------------------------------//

    if (monitor(iteration, r)) break;

  }

  // copy into the solution vector if needed
  if (x0 != &x) x.copy(*x0);

}

} // end namespace callow

#endif // callow_JACOBI_I_HH_

//----------------------------------------------------------------------------//
//              end of file Jacobi.i.hh
//----------------------------------------------------------------------------//
