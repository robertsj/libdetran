//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Vector.i.hh
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  Vector inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_VECTOR_I_HH_
#define callow_VECTOR_I_HH_

#include "utilities/DBC.hh"
#include <cmath>
#include <stdio.h>
#include <iostream>

namespace callow
{

//---------------------------------------------------------------------------//
// ACCESS
//---------------------------------------------------------------------------//

inline const double& Vector::operator[](const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

inline double& Vector::operator[](const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

inline const double& Vector::operator()(const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

inline double& Vector::operator()(const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

inline const double& Vector::value(const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

inline double& Vector::value(const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

//---------------------------------------------------------------------------//
// VECTOR OPERATIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
inline double Vector::norm(const int type)
{
  double val = 0.0;
#ifdef CALLOW_ENABLE_PETSC_OPS
  if (type == L1)
    VecNorm(d_petsc_vector, NORM_1, &val);
  else if (type == L2 || type == L2GRID)
    VecNorm(d_petsc_vector, NORM_2, &val);
  else if (type == LINF)
    VecNorm(d_petsc_vector, NORM_INFINITY, &val);
  else
    THROW("Unsupported norm type");
#else
  if (type == L1 || type == L1GRID)
  {
    for (int i = 0; i < d_size; i++)
      val += std::abs(d_value[i]);
  }
  else if (type == L2 || type == L2GRID)
  {
    for (int i = 0; i < d_size; i++)
      val += d_value[i] * d_value[i];
    val = std::sqrt(val);
  }
  else if (type == LINF)
  {
    for (int i = 0; i < d_size; i++)
      val = std::max(val, std::abs(d_value[i]));
  }
#endif
  // divide by N or sqrt(N) for the grid norms
  if (type == L1GRID) val /= (double)d_size;
  if (type == L2GRID) val /= std::sqrt((double)d_size);
  return val;
}

//---------------------------------------------------------------------------//
inline double Vector::norm_residual(const Vector& x, const int type)
{
  Require(d_size == x.size());
  double val = 0.0;
#ifdef CALLOW_ENABLE_PETSC_OPS
  // create the vector tmp = me - x;
  Vector tmp(*this);
  VecAXPY(tmp.petsc_vector(), -1.0, const_cast<Vector*>(&x)->petsc_vector());
  // and take its norm
  val = tmp.norm(type);
#else
  // basic norms
  if (type == L1)
  {
    for (int i = 0; i < d_size; i++)
      val += std::abs(d_value[i] - x[i]);
  }
  else if (type == L2)
  {
    for (int i = 0; i < d_size; i++)
      val += (d_value[i] - x[i])*(d_value[i] - x[i]);
    val = std::sqrt(val);
  }
  else if (type == LINF)
  {
    for (int i = 0; i < d_size; i++)
      val = std::max(val, std::abs(d_value[i] - x[i]));
  }
  // relative norms
  else if (type == L1REL)
  {
    for (int i = 0; i < d_size; i++)
      val += std::abs((d_value[i] - x[i])/d_value[i]);
  }
  else if (type == L2REL)
  {
    for (int i = 0; i < d_size; i++)
      val += ((d_value[i] - x[i])/d_value[i])*((d_value[i] - x[i])/d_value[i]);
    val = std::sqrt(val);
  }
  else if (type == LINFREL)
  {
    for (int i = 0; i < d_size; i++)
      val = std::max(val, std::abs(d_value[i] - x[i]));
  }
  else
    THROW("Unsupported norm residual type");
#endif
  return val;
}

inline double Vector::norm_residual(SP_vector x, const int type)
{
  Require(x);
  return norm_residual(*x, type);
}

//---------------------------------------------------------------------------//
inline void Vector::set(const double v)
{
#ifdef CALLOW_ENABLE_PETSC_OPS2
  VecSet(d_petsc_vector, v);
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] = v;
#endif
}

//---------------------------------------------------------------------------//
inline void Vector::scale(const double v)
{
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecScale(d_petsc_vector, v);
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] *= v;
#endif
}

//---------------------------------------------------------------------------//
inline double Vector::dot(const Vector& x)
{
  Require(d_size == x.size());
  double val = 0;
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecDot(d_petsc_vector, const_cast<Vector* >(&x)->petsc_vector(), &val);
#else
  for (int i = 0; i < d_size; i++)
    val += d_value[i] * x[i];
#endif
  return val;
}

inline double Vector::dot(SP_vector x)
{
  Require(x);
  return dot(*x);
}

//---------------------------------------------------------------------------//
inline void Vector::add(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecAXPY(d_petsc_vector, 1.0, const_cast<Vector* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] += x[i];
#endif
}

inline void Vector::add(SP_vector x)
{
  Require(x);
  add(*x);
}

//---------------------------------------------------------------------------//
inline void Vector::add_a_times_x(const double a, const Vector& x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS2
  VecAXPY(d_petsc_vector, a, const_cast<Vector* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] += a*x[i];
#endif
}

inline void Vector::add_a_times_x(const double a, SP_vector x)
{
  Require(x);
  add_a_times_x(a, *x);
}

//---------------------------------------------------------------------------//
inline void Vector::subtract(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecAXPY(d_petsc_vector, -1.0, const_cast<Vector* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] -= x[i];
#endif
}

inline void Vector::subtract(SP_vector x)
{
  Require(x);
  subtract(*x);
}

//---------------------------------------------------------------------------//
inline void Vector::multiply(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  Vector tmp(*this);
  VecPointwiseMult(d_petsc_vector, tmp.petsc_vector(), const_cast<Vector* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] *= x[i];
#endif
}

inline void Vector::multiply(SP_vector x)
{
  Require(x);
  multiply(*x);
}

//---------------------------------------------------------------------------//
inline void Vector::divide(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  Vector tmp(*this);
  VecPointwiseDivide(d_petsc_vector, tmp.petsc_vector(), const_cast<Vector* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] /= x[i];
#endif
}

inline void Vector::divide(SP_vector x)
{
  Require(x);
  divide(*x);
}

//---------------------------------------------------------------------------//
inline void Vector::copy(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecCopy(const_cast<Vector* >(&x)->petsc_vector(), d_petsc_vector);
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] = x[i];
#endif
}

inline void Vector::copy(SP_vector x)
{
  Require(x);
  copy(*x);
}

} // end namespace callow

#endif /* callow_VECTOR_I_HH_ */
