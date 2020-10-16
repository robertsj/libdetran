//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MathUtilities.hh
 *  @brief Provides several useful math functions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 *
 *  @todo There is some overlap with callow vector methods.  It might be
 *        a good idea to move all math extras to callow
 */
//----------------------------------------------------------------------------//

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
inline double norm(const vec_dbl &x, std::string flag = "L2")
{
  double v = 0.0;
  if (flag == "L2")
  {
    for (size_t i = 0; i < x.size(); i++)
      v += x[i]*x[i];
    v = std::sqrt(v);
  }
  else if (flag == "L1")
  {
    for (size_t i = 0; i < x.size(); i++)
      v += std::abs(x[i]);
  }
  else if (flag == "Linf")
  {
    v = x[0];
    for (size_t i = 1; i < x.size(); i++)
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
  for (size_t i = 0; i < x.size(); ++i)
    x[i] = x[i] * scale;
}

/// Set a double vector.
inline void vec_set(vec_dbl &x, double value)
{
  for (size_t i = 0; i < x.size(); ++i)
    x[i] = value;
}

/// Norm of the residual of two double vectors.
inline double norm_residual(const vec_dbl &x,
                            const vec_dbl &y,
                            std::string    flag = "L2")
{
  Require(x.size() == y.size());
  double v = 0.0;
  if (flag == "L2")
  {
    // Square root of sum of square of differences.
    for (size_t i = 0; i < x.size(); ++i)
      v += std::pow(x[i]-y[i], 2);
    v = std::sqrt(v);
  }
  else if (flag == "L1")
  {
    // Sum of absolute value of differences.
    for (size_t i = 0; i < x.size(); ++i)
      v += std::abs(x[i]-y[i]);
  }
  else if (flag == "Linf")
  {
    // Max absolute value of differences.
    v = std::abs(x[0] - y[0]);
    for (size_t i = 1; i < x.size(); ++i)
      v = std::max(std::abs(x[i]-y[i]), v);
  }
  else
  {
    THROW("Bad norm residual flag.");
  }
  return v;
}

/// Norm of the relative residual.
inline double norm_relative_residual(const vec_dbl &x,
                                     const vec_dbl &y,
                                     std::string    flag = "L2")
{
  Require(x.size() == y.size());
  double v = 0.0;
  if (flag == "L2")
  {
    // Square root of sum of square of differences.
    for (size_t i = 0; i < x.size(); ++i)
      v += std::pow((x[i]-y[i])/y[i], 2);
    v = std::sqrt(v);
  }
  else if (flag == "L1")
  {
    // Sum of absolute value of differences.
    for (size_t i = 0; i < x.size(); ++i)
      v += std::abs((x[i]-y[i])/y[i]);
  }
  else if (flag == "Linf")
  {
    // Max absolute value of differences.
    v = std::abs(x[0] - y[0]);
    for (size_t i = 1; i < x.size(); ++i)
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
  for (size_t i = 0; i < v.size(); ++i)
    if (v[i] < m) m = v[i];
  return m;
}

/// Return the maximum of a vector.
template <class T>
inline T vec_max(const std::vector<T> &v)
{
  T m = std::numeric_limits<T>::min();
  for (size_t i = 0; i < v.size(); ++i)
    if (v[i] > m) m = v[i];
  return m;
}

/// Return the sum of a vector.
template <class T>
inline T vec_sum(const std::vector<T> &v)
{
  T s = 0;
  for (size_t i = 0; i < v.size(); ++i)
    s += v[i];
  return s;
}

/// Add a constant to a vector
template <class T>
inline void vec_plus_a(std::vector<T> &v, const T a)
{
  for (size_t i = 0; i < v.size(); ++i)
    v[i] += a;
}

/// Return the unique elements of a vector
template <class T>
inline std::vector<T> vec_unique(const std::vector<T> &v)
{
  std::vector<T> u = v;                   // 1 1 2 3 5 6 6 1
  std::sort(u.begin(), u.end());          // 1 1 1 2 3 5 6 6
  typename std::vector<T>::iterator it;
  it = std::unique(u.begin(), u.end());   // 1 2 3 5 6 ? ? ?
  u.resize(std::distance(u.begin(), it)); // 1 2 3 5 6
  return u;
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
  std::vector<double> v(n, a + 0.5*dx);
  for (int i = 1; i < n; ++i)
    v[i] = v[i-1] + dx;
  return v;
}

/**
 *  @brief Create a range vector
 *
 *  Example: range(0, 4)       = {0, 1, 2, 3}
 *           range(4, 0)       = {4, 3, 2, 1}
 *           range(0, 4, true) = {0, 1, 2, 3, 4}
 *
 *  @param  a    start of range
 *  @param  b    end of range
 */
template <class T>
inline std::vector<T> range(const int a, const int b, bool inclusive = false)
{
  int size = std::abs(a-b);
  if (inclusive) ++size;
  std::vector<T> v(size, 0);
  typename std::vector<T>::iterator it = v.begin();
  int del = a < b ? 1 : -1;
  int i = a;
  for (; it != v.end(); ++it, i += del)
  {
    *it = i;
  }
  return v;
}


} // namespace detran_utilities

#endif /* detran_utilities_MATHUTILITIES_HH_ */

//----------------------------------------------------------------------------//
//              end of MathUtilities.hh
//----------------------------------------------------------------------------//
