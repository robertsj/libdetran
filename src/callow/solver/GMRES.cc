//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GMRES.cc
 * \author robertsj
 * \date   Sep 14, 2012
 * \brief  GMRES class definition.
 */
//---------------------------------------------------------------------------//

#include "GMRES.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

GMRES::GMRES(const double  atol,
                const double  rtol,
                const int     maxit,
                const int     restart)
  : LinearSolver(atol, rtol, maxit, "solver_gmres")
  , d_restart(restart)
  , d_c(restart+1, 0.0)
  , d_s(restart+1, 0.0)
  , d_reorthog(1)
{
  Insist(d_restart > 2, "Need a restart of > 2");
  d_H = new double*[(restart + 1)];
  for (int i = 0; i <= d_restart; i++)
  {
    d_H[i] = new double[restart];
    for (int j = 0; j < d_restart; j++)
      d_H[i][j] = 0.0;
  }
//  d_c.resize(restart + 1);
//  d_s.resize(restart + 1);
}


GMRES::~GMRES()
{
  for (int i = 0; i <= d_restart; i++)
  {
    delete [] d_H[i];
  }
  delete [] d_H;
}

} // end namespace callow



