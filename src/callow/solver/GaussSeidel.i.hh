//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussSeidel.i.hh
 *  @author robertsj
 *  @date   Sep 14, 2012
 *  @brief  GaussSeidel.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_GAUSSSEIDEL_I_HH_
#define callow_GAUSSSEIDEL_I_HH_

#include "callow/matrix/Matrix.hh"
#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

inline void GaussSeidel::solve_impl(const Vector &b, Vector &x)
{

  Insist(dynamic_cast< Matrix* >(d_A.bp()),
    "Need an explicit matrix for use with GaussSeidel iteration");
  Matrix::SP_matrix A = d_A;

  // temporary storage and pointers for swapping
  Vector temp(x.size(), 0.0);
  Vector temp2(x.size(), 0.0);
  Vector* x0 = &x;
  Vector* x1 = &temp;
  Vector* swap;

  // compute initial residual Ax - b and its norm
  A->multiply((*x0), (*x1));
  double r = x1->norm_residual(b, L2);
  if (monitor_init(r)) 
  {
    //return;
  }
  // perform iterations
  for (int iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- -inv(D+L)*U*X0 + inv(D+L)*b
    //---------------------------------------------------//

    double* a = A->values();
    for (int i = 0; i < A->number_rows(); i++)
    {
      double v = 0;
      int p = A->start(i);
      int d = A->diagonal(i);
      // lower triangle -- we have updated unknowns
      for (; p < d; ++p)
        v += a[p] * (*x1)[A->column(p)];
      ++p; // skip diagonal
      // upper triangle -- we do not have updated unknowns
      for (; p < A->end(i); ++p)
        v += a[p] * (*x0)[A->column(p)];
      (*x1)[i] = (b[i] - v) / a[d];
    }
    a = 0; // nullify pointer

    //---------------------------------------------------//
    // relax
    //---------------------------------------------------//

    if (d_omega != 1.0)
    {
      // x1 <-- (1-w)*x0 + w*x1
      x1->scale(d_omega);
      x1->add_a_times_x((1.0-d_omega), *x0);
    }

    //---------------------------------------------------//
    // compute residual norm
    //---------------------------------------------------//

    if (d_successive_norm)
    {
      // compare x1 with x0
      r = x1->norm_residual(*x0, d_norm_type);
    }
    else
    {
      // Compute the residual and put it into x0
      A->multiply(*x1, temp2);
      r = temp2.norm_residual(b, d_norm_type);
    }

    //---------------------------------------------------//
    // swap pointers
    //---------------------------------------------------//
    swap = x0;
    x0   = x1;
    x1   = swap;

    //---------------------------------------------------//
    // check convergence
    //---------------------------------------------------//

    if (monitor(iteration, r)) break;

  }

  // copy into the solution vector if needed
  if (x0 != &x) x.copy(*x0);

}

} // end namespace callow


#endif /* callow_GAUSSSEIDEL_I_HH_ */
