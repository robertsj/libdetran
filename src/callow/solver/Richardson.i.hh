//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Richardson.i.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Richardson inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_RICHARDSON_I_HH_
#define callow_RICHARDSON_I_HH_

#include <cmath>
#include <cstdio>

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
Richardson<T>::Richardson(const double  atol,
                          const double  rtol,
                          const int     maxit,
                          const double  omega)
  : LinearSolver<T>(atol, rtol, maxit, "solver_richardson")
  , d_omega(omega)
{

}

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//


template <class T>
inline void Richardson<T>::solve_impl(const Vector<T> &b, Vector<T> &x)
{

  typedef Vector<T> Vec;

  // scale rhs by relaxation factor
  Vec B(b.size(), 0.0);
  B.add(b);
  B.scale(d_omega);

  // temporary storage and pointers for swapping
  Vec temp(x.size(), 0.0);
  Vec* x0 = &x;
  Vec* x1 = &temp;
  Vec* swap;

  // compute initial residual w(Ax - b) and its norm
  d_A->multiply((*x0), (*x1));
  x1->scale(d_omega);
  T r = x1->norm_residual(B, Vec::L2);
  if (monitor_init(r)) return;

  // perform iterations
  for (int iteration = 1; iteration <= d_maximum_iterations; ++iteration)
  {

    //---------------------------------------------------//
    // compute X1 <-- (I - w*A) * X0 + w*b
    //---------------------------------------------------//

    // X1 <-- A * X0
    d_A->multiply((*x0), (*x1));
    // X1 <-- w * X1 = w * A * X0
    x1->scale(d_omega);
    // X1 <-- X1 - X0 =  (A - I) * X0
    x1->subtract(*x0);
    // X1 <-- X1 - b = (A - I) * X0 - b
    x1->subtract(B);
    // X1 <-- -X1 = (I - A) * X0 + b
    x1->scale(-1);

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

  return;
}

} // end namespace callow

#endif /* callow_RICHARDSON_I_HH_ */
