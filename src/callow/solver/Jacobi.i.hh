//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Jacobi.i.hh
 * \brief  Jacobi inline member definitions
 * \author Jeremy Roberts
 * \date   Sep 13, 2012
 */
//---------------------------------------------------------------------------//

#ifndef callow_JACOBI_I_HH_
#define callow_JACOBI_I_HH_

#include "callow/matrix/Matrix.hh"
#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Jacobi<T>::Jacobi(const double  atol,
                  const double  rtol,
                  const int     maxit)
  : LinearSolver<T>(atol, rtol, maxit, "solver_jacobi")
{
  /* ... */
}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

template <class T>
inline void Jacobi<T>::solve_impl(const Vector<T> &b, Vector<T> &x)
{

  Insist(dynamic_cast< Matrix<T>* >(d_A.bp()),
    "Need an explicit matrix for use with Jacobi iteration");
  typename Matrix<T>::SP_matrix A = d_A;

  typedef Vector<T> Vec;

  // temporary storage and pointers for swapping
  Vec temp(x.size(), 0.0);
  Vec* x0 = &x;
  Vec* x1 = &temp;
  Vec* swap;

  // iteration count
  int &iteration = d_number_iterations;

  // compute initial residual Ax - b and its norm
  A->multiply((*x0), (*x1));
  T r = x1->norm_residual(b, L2);
  if (monitor_init(r))
  {
    //return;
  }

  // perform iterations
  for (int iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- -inv(D)*(L+U)*X0 + inv(D)*b
    //---------------------------------------------------//

    T* a = A->values();
    for (int i = 0; i < A->number_rows(); i++)
    {
      T v = 0;
      int p = A->start(i);
      int d = A->diagonal(i);
      // L * X0
      for (; p < d; ++p)
        v += a[p] * (*x0)[A->column(p)];
      ++p; // skip diagonal
      // U * X0
      for (; p < A->end(i); ++p)
        v += a[p] * (*x0)[A->column(p)];
      (*x1)[i] = (b[i] - v) / a[d];
    }
    a = 0; // nullify pointer

    //---------------------------------------------------//
    // compute residual norm
    //---------------------------------------------------//

    r = x1->norm_residual(*x0, L2);

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

#endif // callow_JACOBI_I_HH_

//---------------------------------------------------------------------------//
//              end of file Jacobi.i.hh
//---------------------------------------------------------------------------//
