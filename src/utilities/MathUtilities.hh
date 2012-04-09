//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MathUtilities.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  MathUtilities class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef MATHUTILITIES_HH_
#define MATHUTILITIES_HH_

// Utilities
#include "DBC.hh"
#include "Definitions.hh"

// System
#include <cmath>

namespace detran
{

/// L2 norm of the residual of two vectors.
inline double norm(vec_dbl &x)
{
  double v = 0.0;
  for (int i = 0; i < x.size(); i++)
    v += pow(x[i]*x[i], 2);
  return std::sqrt(v);
}

inline double norm_residual(vec_dbl &x, vec_dbl &y)
{
  Require(x.size() == y.size());
  double v = 0.0;
  for (int i = 0; i < x.size(); i++)
    v += pow(x[i]-y[i], 2);
  return std::sqrt(v);
}

inline double norm_relative_residual(vec_dbl &x, vec_dbl &y)
{
  Require(x.size() == y.size());
  double v = 0.0;
  for (int i = 0; i < x.size(); i++)
    v += pow((x[i]-y[i])/y[i], 2);
  return std::sqrt(v);
}

} // namespace

#endif /* MATHUTILITIES_HH_ */

//---------------------------------------------------------------------------//
//              end of MathUtilities.hh
//---------------------------------------------------------------------------//
