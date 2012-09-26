//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GMRES.i.hh
 *  @author robertsj
 *  @date   Sep 14, 2012
 *  @brief  GMRES inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_GMRES_I_HH_
#define callow_GMRES_I_HH_

#include "matrix/Matrix.hh"
#include <cmath>
#include <cstdio>
#include <iostream>
using std::cout;
using std::endl;

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

template <class T>
GMRES<T>::GMRES(const double  atol,
                const double  rtol,
                const int     maxit,
                const int     restart)
  : LinearSolver<T>(atol, rtol, maxit, "solver_gmres")
  , d_restart(restart)
  , d_c(restart+1, 0.0)
  , d_s(restart+1, 0.0)
  , d_reorthog(1)
{
  Insist(d_restart > 2, "Need a restart of > 2");
  d_H = new T*[(restart + 1)];
  for (int i = 0; i <= d_restart; i++)
  {
    d_H[i] = new T[restart];
    for (int j = 0; j < d_restart; j++)
      d_H[i][j] = 0.0;
  }
//  d_c.resize(restart + 1);
//  d_s.resize(restart + 1);
}

template <class T>
GMRES<T>::~GMRES()
{
  for (int i = 0; i <= d_restart; i++)
  {
    delete [] d_H[i];
  }
  delete [] d_H;
}
//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

template <class T>
inline void GMRES<T>::solve_impl(const Vector<T> &b, Vector<T> &x)
{

  typedef Vector<T> Vec_T;

  int restart = d_restart;
  if (restart >= d_A->number_rows()) restart = d_A->number_rows();

  // krylov basis
  std::vector<Vec_T>  v(d_restart + 1, Vec_T(x.size(), 0.0));

  // residual
  Vec_T r(x.size(), 0.0);
  Vec_T t(x.size(), 0.0);

  // vector such that x = V*y
  Vec_T y(d_restart, 0.0);

  // vector such that g(1:k) = R*y --> x = V*inv(R)*g and |g(k+1)| is the residual
  Vec_T g(d_restart + 1, 0.0);

  // initialize c and s
  d_c.set(0.0);
  d_s.set(0.0);

  //-------------------------------------------------------------------------//
  // outer iterations
  //-------------------------------------------------------------------------//

  int iteration = 0; // outer iteration
  bool done = false;
  while (!done and iteration < d_maximum_iterations)
  {

    // clear krylov subspace
    for (int i = 0; i < d_restart; i++) v[i].set(0.0);
    g.set(0.0);

    // compute residual
    //   apply operator
    d_A->multiply(x, r);
    r.subtract(b);
    r.scale(-1);
    //   apply left preconditioner
    if (d_P and d_pc_side == Base::LEFT)
    {
      t.copy(r);
      d_P->apply(t, r);
    }
    //   compute norm of residual
    T rho = r.norm(L2);

    // otherwise, we must be restarting

    // check initial outer residual.  if it's small enough, we started
    // with a solved system.
    if (iteration == 0)
    {
      //cout << "initial..." << endl;
      if (monitor_init(rho)) return;
    }
    else
    {
      //if (d_monitor_output) cout << "restarting..." << endl;
    }

    // initial krylov vector
    v[0].copy(r);
    v[0].scale(1.0 / rho);
    g[0] = rho;

    // inner iterations (of size restart)
    int k = 0;
    int lala = 0;
    for (; k < d_restart; ++k)
    {
      ++iteration;
      // check iteration count
      if (iteration >= d_maximum_iterations-1)
      {
        done = true;
        break;
      }

      //---------------------------------------------------------------------//
      // compute v(k+1) <-- inv(P_L)*A*inv(P_R) * v(k)
      //---------------------------------------------------------------------//

      // right preconditioner
      if (d_P and d_pc_side == Base::RIGHT)
      {
        t.copy(v[k]);
        d_P->apply(t, v[k]);
      }
      // apply A
      d_A->multiply(v[k], v[k+1]);
      // left preconditioner
      if (d_P and d_pc_side == Base::LEFT)
      {
        t.copy(v[k+1]);
        d_P->apply(t, v[k+1]);
      }

      //---------------------------------------------------------------------//
      // use modified gram-schmidt to orthogonalize v(k+1)
      //---------------------------------------------------------------------//

      T norm_Av = v[k+1].norm();
      for (int j = 0; j <= k; ++j)
      {
        d_H[j][k] = v[k+1].dot(v[j]);
        v[k+1].add_a_times_x(-d_H[j][k], v[j]);
      }
      d_H[k+1][k] = v[k+1].norm(L2);
      T norm_Av_2 = d_H[k+1][k];

      //---------------------------------------------------------------------//
      // optional reorthogonalization
      //---------------------------------------------------------------------//

      if ( (d_reorthog == 1 and norm_Av + 0.001 * norm_Av_2 == norm_Av) or
           (d_reorthog == 2) )
      {
        // summarized from kelley:
        //  if the new vector (i.e. v[k+1]) is very small relative to
        //  A*v[k], then information might be lost so reorthogonalize.  the
        //  delta of 0.001 is what kelley uses in his test code.

        cout << " reorthog ... " << endl;
        for (int j = 0; j < k; ++j)
        {
          T hr = v[j].dot(v[k+1]);
          d_H[j][k] += hr;
          v[k+1].add_a_times_x(-hr, v[j]);
        }
        d_H[k+1][k] = v[k+1].norm();
      }

      //---------------------------------------------------------------------//
      // watch for happy breakdown: if H[k+1][k] == 0, we've solved Ax=b
      //---------------------------------------------------------------------//

      if (d_H[k+1][k] != 0.0)
      {
        v[k+1].scale(1.0/d_H[k+1][k]);
      }
      else
      {
        std::printf("happy breakdown for k = %5i (iteration = %5i) \n",
                    k, iteration);
      }

      //---------------------------------------------------------------------//
      // apply givens rotations to triangularize H on-the-fly (it's neat!)
      //---------------------------------------------------------------------//

      if (k > 0) apply_givens(k);
      double nu = std::sqrt(d_H[k][k]*d_H[k][k] + d_H[k+1][k]*d_H[k+1][k]);
      d_c[k] =  d_H[k  ][k] / nu;
      d_s[k] = -d_H[k+1][k] / nu;
      d_H[k  ][k] = d_c[k] * d_H[k][k] - d_s[k]*d_H[k+1][k];
      d_H[k+1][k] = 0.0;
      T g_0 = d_c[k]*g[k] - d_s[k]*g[k+1];
      T g_1 = d_s[k]*g[k] + d_c[k]*g[k+1];
      g[k  ] = g_0;
      g[k+1] = g_1;

      //---------------------------------------------------------------------//
      // monitor the residual and break if done
      //---------------------------------------------------------------------//

      rho = std::abs(g_1);
      if (monitor(iteration, rho))
      {
        ++k;
//        printf("gmres(%3i) terminated at outer iteration %5i ",
//               d_restart, iteration/d_restart);
//        printf("(inner iteration %3i) to a solution with residual: %12.8e \n",
//               k, rho);
        done = true;
        break;
      }

    } // end inners

    //---------------------------------------------------------------------//
    // update the solution
    //---------------------------------------------------------------------//

    compute_y(y, g, k);

    // update x = x0 + v[0]*y[0] + ...
    for (int i = 0; i < k; ++i)
    {
      x.add_a_times_x(y[i], v[i]);
    }
    // \todo this assumes x_0 = 0; otherwise, need x = x_0 + inv(P)*(V*y)
    if (d_P and d_pc_side == Base::RIGHT)
    {
      t.copy(x);
      d_P->apply(t, x);
    }


  } // end outers

  d_A->multiply(x, r);
  r.subtract(b);
  r.scale(-1.0);
  T resid = r.norm(L2);
  //printf(" final residual norm is actually: %12.8e \n", resid);

  return;

}

template <class T>
inline void GMRES<T>::apply_givens(const int k)
{
  Require(k < d_restart);
  for (int i = 0; i < k; ++i)
  {
    T g_0 = d_c[i]*d_H[i][k] - d_s[i]*d_H[i+1][k];
    T g_1 = d_s[i]*d_H[i][k] + d_c[i]*d_H[i+1][k];
    d_H[i][k]   = g_0;
    d_H[i+1][k] = g_1;
  }
}

template <class T>
inline void GMRES<T>::compute_y(Vector<T> &y, const Vector<T> &g, const int k)
{
  // H is [m+1][m], though we may have only need for k*k
  // solves H[0:k][0:k]*y[0:k] = g[0:k]
  //  but only for H_ij for j>=i, i.e. upper triangle
  Require(k <= d_restart);
  Require(y.size() >= k);
  Require(g.size() >= k);

  for (int i = k - 1; i >= 0; --i)
  {
    y[i] = g[i];
    for (int j = i + 1; j < k; ++j)
    {
      y[i] -= d_H[i][j] * y[j];
    }
    Assert(d_H[i][i] != 0.0);
    y[i] /= d_H[i][i];
  }

}

template <class T>
inline void GMRES<T>::initialize_H()
{
  for (int i = 0; i <= d_restart; i++)
    for (int j = 0; j < d_restart; j++)
      d_H[i][j] = 0.0;
}

} // end namespace callow

#endif /* callow_GMRES_I_HH_ */
