//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MathUtilities.hh
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  Provides several useful math functions
 */
//---------------------------------------------------------------------------//

#ifndef detran_utilities_MATHUTILITIES_HH_
#define detran_utilities_MATHUTILITIES_HH_

#include "DBC.hh"
#include "Definitions.hh"

#include <algorithm>
#include <cmath>
#include <string>
#include <limits>

namespace detran_utilities
{

/// Norm of a vector
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


/// Norm of the residual of two double vectors.
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

/// Norm of the relative residual.
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

/// Return the minimum of a vector.
template <class T>
inline T vec_min(const std::vector<T> &v)
{
  T m = std::numeric_limits<T>::max();
  for (int i = 0; i < v.size(); ++i)
    if (v[i] < m) m = v[i];
  return m;
}

/// Return the maximum of a vector.
template <class T>
inline T vec_max(const std::vector<T> &v)
{
  T m = std::numeric_limits<T>::min();
  for (int i = 0; i < v.size(); ++i)
    if (v[i] > m) m = v[i];
  return m;
}

/// Linear mesh between a and b.  Gives the edges.
inline std::vector<double> linspace(double a, double b, int n = 10)
{
  Require(a < b);
  std::vector<double> v(n, a);
  double dx;
  if (n <= 1)
    dx = b-a;
  else
    dx = (b-a)/double(n-1);
  for (int i = 1; i < n; ++i)
    v[i] = v[i-1] + dx;
  return v;
}

/// Linear mesh between a and b. Gives the mesh centers.
inline std::vector<double> linspace_center(double a, double b, int n = 10)
{
  Require(a < b);
  double dx;
  if (n <= 1)
    dx = b-a;
  else
    dx = (b-a)/double(n);
  std::vector<double> v(n, 0.5*dx);
  for (int i = 1; i < n; ++i)
    v[i] = v[i-1] + dx;
  return v;
}

} // namespace detran_utilities

#endif /* detran_utilities_MATHUTILITIES_HH_ */

//---------------------------------------------------------------------------//
//              end of MathUtilities.hh
//---------------------------------------------------------------------------//
