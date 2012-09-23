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
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
GaussSeidel<T>::GaussSeidel(const double  atol,
                  const double  rtol,
                  const int     maxit)
  : LinearSolver<T>(atol, rtol, maxit, "solver_gauss_seidel")
{
  /* ... */
}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

template <class T>
inline void GaussSeidel<T>::solve_impl(const Vector<T> &b, Vector<T> &x)
{

  Insist(dynamic_cast< Matrix<T>* >(d_A.bp()),
    "Need an explicit matrix for use with GaussSeidel iteration");
  typename Matrix<T>::SP_matrix A = d_A;

  typedef Vector<T> Vec;

  // temporary storage and pointers for swapping
  Vec temp(x.size(), 0.0);
  Vec* x0 = &x;
  Vec* x1 = &temp;
  Vec* swap;

  // compute initial residual Ax - b and its norm
  A->multiply((*x0), (*x1));
  T r = x1->norm_residual(b, Vec::L2);
  if (monitor_init(r)) return;

  // perform iterations
  for (int iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- -inv(D+L)*U*X0 + inv(D+L)*b
    //---------------------------------------------------//

    T* a = A->values();
    for (int i = 0; i < A->number_rows(); i++)
    {
      T v = 0;
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
    // compute residual norm
    //---------------------------------------------------//

    r = x1->norm_residual(*x0, d_norm_type);

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
