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
#include <algorithm>
#include <cmath>

namespace detran
{

/// L2 norm of the residual of two vectors.
inline double norm(vec_dbl &x)
{
  double v = 0.0;
  for (int i = 0; i < x.size(); i++)
    v += std::pow(x[i]*x[i], 2);
  return std::sqrt(v);
}

/*!
 *  \brief Norm of the residual of two double vectors.
 *  \param  flag    Default is false for L2.  True for L-infinity.
 */
inline double norm_residual(vec_dbl &x, vec_dbl &y, bool flag = false)
{
  Require(x.size() == y.size());
  double v = 0.0;
  // L2 norm
  if (!flag)
  {
    for (int i = 0; i < x.size(); i++)
      v += std::pow(x[i]-y[i], 2);
    return std::sqrt(v);
  }
  // L-infinity norm
  else
  {
    v = std::abs(x[0] - y[0]);
    for (int i = 1; i < x.size(); i++)
      v = std::max(std::abs(x[i]-y[i]), v);
    return v;
  }
}

inline double norm_relative_residual(vec_dbl &x, vec_dbl &y)
{
  Require(x.size() == y.size());
  double v = 0.0;
  for (int i = 0; i < x.size(); i++)
    v += std::pow((x[i]-y[i])/y[i], 2);
  return std::sqrt(v);
}

} // namespace

#endif /* MATHUTILITIES_HH_ */

//---------------------------------------------------------------------------//
//              end of MathUtilities.hh
//---------------------------------------------------------------------------//
