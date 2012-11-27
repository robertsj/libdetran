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
// IMPLEMENTATION
//---------------------------------------------------------------------------//

inline void GMRES::solve_impl(const Vector &b, Vector &x0)
{
  int restart = d_restart;
  if (restart >= d_A->number_rows()) restart = d_A->number_rows();

  // Unknowns.  If x0 is nonzero, we need to separate it out.
  Vector x(x0);

  // krylov basis
  std::vector<Vector>  v(d_restart + 1, Vector(x.size(), 0.0));

  // residual
  Vector r(x.size(), 0.0);
  Vector t(x.size(), 0.0);

  // vector such that x = V*y
  Vector y(d_restart, 0.0);

  // vector such that g(1:k) = R*y --> x = V*inv(R)*g and |g(k+1)| is the residual
  Vector g(d_restart + 1, 0.0);

  // initialize c and s
  d_c.set(0.0);
  d_s.set(0.0);

  //-------------------------------------------------------------------------//
  // outer iterations
  //-------------------------------------------------------------------------//

  int iteration = 0;  // total iterations (i.e. applications of A)
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
    double rho = r.norm(L2);

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

      // apply right preconditioner and operator
      if (d_P and d_pc_side == Base::RIGHT)
      {
        d_P->apply(v[k], v[k + 1]);
        t.copy(v[k + 1]);
        d_A->multiply(t, v[k + 1]);
      }
      else
      {
        // save on a copy
        d_A->multiply(v[k], v[k + 1]);
      }
      // apply left preconditioner
      if (d_P and d_pc_side == Base::LEFT)
      {
        t.copy(v[k + 1]);
        d_P->apply(t, v[k + 1]);
      }

      //---------------------------------------------------------------------//
      // use modified gram-schmidt to orthogonalize v(k+1)
      //---------------------------------------------------------------------//

      double norm_Av = v[k+1].norm();
      for (int j = 0; j <= k; ++j)
      {
        d_H[j][k] = v[k+1].dot(v[j]);
        v[k+1].add_a_times_x(-d_H[j][k], v[j]);
      }
      d_H[k+1][k] = v[k+1].norm(L2);
      double norm_Av_2 = d_H[k+1][k];

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

        if (d_monitor_level > 1) cout << " reorthog ... " << endl;
        for (int j = 0; j < k; ++j)
        {
          double hr = v[j].dot(v[k+1]);
          d_H[j][k] += hr;
          v[k+1].add_a_times_x(-hr, v[j]);
        }
        d_H[k+1][k] = v[k+1].norm();
      }

      //---------------------------------------------------------------------//
      // watch for happy breakdown: if H[k+1][k] == 0, we've solved Ax=b
      //---------------------------------------------------------------------//
      bool happy = false;
      if (d_H[k+1][k] != 0.0)
      {
        v[k+1].scale(1.0/d_H[k+1][k]);
      }
      else
      {
        happy = true;
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
      double g_0 = d_c[k]*g[k] - d_s[k]*g[k+1];
      double g_1 = d_s[k]*g[k] + d_c[k]*g[k+1];
      g[k  ] = g_0;
      g[k+1] = g_1;
      if (happy) cout << " g_0 = " << g_0 << " g_1 = " << g_1 << endl;

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

    // reset the solution
    x.set(0.0);
    // update x = v[0]*y[0] + ...
    for (int i = 0; i < k; ++i)
    {
      x.add_a_times_x(y[i], v[i]);
    }
    if (d_P and d_pc_side == Base::RIGHT)
    {
      t.copy(x);
      d_P->apply(t, x);
    }
    // add the initial guess
    x.add(x0);

    // This iterate becomes the new initial guess
    x0.copy(x);

  } // end outers

  // copy solution to outgoing vector
  x0.copy(x);

  return;
}

inline void GMRES::apply_givens(const int k)
{
  Require(k < d_restart);
  for (int i = 0; i < k; ++i)
  {
    double g_0 = d_c[i]*d_H[i][k] - d_s[i]*d_H[i+1][k];
    double g_1 = d_s[i]*d_H[i][k] + d_c[i]*d_H[i+1][k];
    d_H[i][k]   = g_0;
    d_H[i+1][k] = g_1;
  }
}

inline void GMRES::compute_y(Vector &y, const Vector &g, const int k)
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

inline void GMRES::initialize_H()
{
  for (int i = 0; i <= d_restart; i++)
    for (int j = 0; j < d_restart; j++)
      d_H[i][j] = 0.0;
}

} // end namespace callow

#endif /* callow_GMRES_I_HH_ */
