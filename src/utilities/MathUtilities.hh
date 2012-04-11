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
#include <string>

namespace detran
{

/*!
 *  \brief Norm of a vector.
 *  \param  flag    L2 (Default), L1, or Linf
 */
inline double norm(vec_dbl &x, std::string flag = "L2")
{
  double v = 0.0;
  if (flag == "L2")
  {
    for (int i = 0; i < x.size(); i++)
      v += x[i]*x[i];
    v = std::sqrt(v);
  }
  else if (flag == "L1")
  {
    for (int i = 0; i < x.size(); i++)
      v += std::abs(x[i]);
  }
  else if (flag == "Linf")
  {
    v = x[0];
    for (int i = 1; i < x.size(); i++)
      v = std::max(std::abs(x[i]), v);
  }
  else
  {
    THROW("Bad norm flag.");
  }
  return v;
}

/// Scale a double vector.
inline void vec_scale(vec_dbl &x, double scale)
{
  for (int i = 0; i < x.size(); i++)
    x[i] = x[i] * scale;
}


/*!
 *  \brief Norm of the residual of two double vectors.
 *  \param  flag    Default is false for L2.  True for L-infinity.
 */
inline double norm_residual(vec_dbl &x, vec_dbl &y, std::string flag = "L2")
{
  Require(x.size() == y.size());
  double v = 0.0;
  if (flag == "L2")
  {
    // Square root of sum of square of differences.
    for (int i = 0; i < x.size(); i++)
      v += std::pow(x[i]-y[i], 2);
    v = std::sqrt(v);
  }
  else if (flag == "L1")
  {
    // Sum of absolute value of differences.
    for (int i = 0; i < x.size(); i++)
      v += std::abs(x[i]-y[i]);
  }
  else if (flag == "Linf")
  {
    // Max absolute value of differences.
    v = std::abs(x[0] - y[0]);
    for (int i = 1; i < x.size(); i++)
      v = std::max(std::abs(x[i]-y[i]), v);
  }
  else
  {
    THROW("Bad norm residual flag.");
  }
  return v;
}

inline double norm_relative_residual(vec_dbl &x, vec_dbl &y, std::string flag = "L2")
{
  Require(x.size() == y.size());
  double v = 0.0;
  if (flag == "L2")
  {
    // Square root of sum of square of differences.
    for (int i = 0; i < x.size(); i++)
      v += std::pow((x[i]-y[i])/y[i], 2);
    v = std::sqrt(v);
  }
  else if (flag == "L1")
  {
    // Sum of absolute value of differences.
    for (int i = 0; i < x.size(); i++)
      v += std::abs((x[i]-y[i])/y[i]);
  }
  else if (flag == "Linf")
  {
    // Max absolute value of differences.
    v = std::abs(x[0] - y[0]);
    for (int i = 1; i < x.size(); i++)
      v = std::max(std::abs((x[i]-y[i])/y[i]), v);
  }
  else
  {
    THROW("Bad norm relative residual flag.");
  }
  return v;
}

} // namespace

#endif /* MATHUTILITIES_HH_ */

//---------------------------------------------------------------------------//
//              end of MathUtilities.hh
//---------------------------------------------------------------------------//
