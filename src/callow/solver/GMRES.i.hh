//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  GMRES.i.hh
 *  @brief GMRES inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

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

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

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
