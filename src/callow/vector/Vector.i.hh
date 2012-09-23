//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Vector.i.hh
 * \author robertsj
 * \date   Sep 13, 2012
 * \brief  Vector inline member definitions
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
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <class T>
Vector<T>::Vector()
  : d_size(0)
  , d_temporary(false)
{
  /* ... */
}

//---------------------------------------------------------------------------//
template <class T>
Vector<T>::Vector(const int n, T v)
  : d_size(n)
  , d_temporary(false)
{
  // Preconditions
  Require(d_size > 0);

  // Size the vector
#ifdef CALLOW_ENABLE_PETSC
  // Create the vector and initialize values
  PetscErrorCode ierr;
  ierr = VecCreateSeq(PETSC_COMM_SELF, n, &d_petsc_vector);
  ierr = VecSet(d_petsc_vector, v);
  // Grab the underlying storage
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  d_value = new T[n];
  set(v);
#endif
  // Postconditions
  Ensure(d_value);
}

//---------------------------------------------------------------------------//
template <class T>
Vector<T>::Vector(const Vector &x)
  : d_size(x.size())
  , d_temporary(false)
{
  // nothing to do if x has no elements
  if (!d_size) return;
#ifdef CALLOW_ENABLE_PETSC
  // create the vector and initialize values
  PetscErrorCode ierr;
  Vec x_v = const_cast<Vector<T>* >(&x)->petsc_vector();
  ierr = VecDuplicate(x_v, &d_petsc_vector);
  ierr = VecCopy(x_v, d_petsc_vector);
  // Grab the underlying storage
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  d_value = new T[d_size];
  set(0.0);
  add(x);
#endif
}

//---------------------------------------------------------------------------//
template <class T>
Vector<T>::Vector(Vector &x)
  : d_size(x.size())
  , d_temporary(false)
{
  // nothing to do if x has no elements
  if (!d_size) return;
#ifdef CALLOW_ENABLE_PETSC
  // create the vector and initialize values
  PetscErrorCode ierr;
  Vec x_v = x.petsc_vector();
  ierr = VecDuplicate(x_v, &d_petsc_vector);
  ierr = VecCopy(x_v, d_petsc_vector);
  // Grab the underlying storage
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  d_value = new T[d_size];
  set(0.0);
  add(x);
#endif
}

//---------------------------------------------------------------------------//
template <class T>
Vector<T>::Vector(std::vector<T> &x)
  : d_size(0)
  , d_temporary(true)
{
  d_value = &x[0];
  d_size  = x.size();
}

//---------------------------------------------------------------------------//
#ifdef CALLOW_ENABLE_PETSC
template <class T>
Vector<T>::Vector(Vec pv)
  : d_size(0)
  , d_temporary(true)
{
  PetscErrorCode ierr;
  ierr = VecGetArray(pv, &d_value);
  ierr = VecRestoreArray(pv, PETSC_NULL);
  ierr = VecGetSize(pv, &d_size);
  d_petsc_vector = pv;
  Ensure(!ierr);
}
#endif

//---------------------------------------------------------------------------//
template <class T>
Vector<T>::~Vector()
{
  if (!d_size or d_temporary) return;
#ifdef CALLOW_ENABLE_PETSC
  d_value = 0;
  VecDestroy(&d_petsc_vector);
#else
  delete [] d_value;
#endif
}

template <class T>
void Vector<T>::resize(const int n, const T v)
{
  Insist(!d_temporary, "Cannot resize a temporary Vector!");
#ifdef CALLOW_ENABLE_PETSC
  // destroy the vector if it's built
  if (d_size)
  {
    d_value = 0;
    VecDestroy(&d_petsc_vector);
  }
  d_size = n;
  // create a new vector
  VecCreateSeq(PETSC_COMM_SELF, n, &d_petsc_vector);
  // grab the underlying storage
  PetscErrorCode ierr;
  ierr = VecGetArray(d_petsc_vector, &d_value);
  ierr = VecRestoreArray(d_petsc_vector, PETSC_NULL);
#else
  if (d_size) delete [] d_value;
  d_size = n;
  if (d_size) d_value = new T[d_size];
  set(v);
#endif
}


//---------------------------------------------------------------------------//
// ACCESS
//---------------------------------------------------------------------------//

template <class T>
inline const T& Vector<T>::operator[](const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

template <class T>
inline T& Vector<T>::operator[](const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

template <class T>
inline const T& Vector<T>::operator()(const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

template <class T>
inline T& Vector<T>::operator()(const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

template <class T>
inline const T& Vector<T>::value(const int i) const
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

template <class T>
inline T& Vector<T>::value(const int i)
{
  Require(i >= 0);
  Require(i < d_size);
  return d_value[i];
}

//---------------------------------------------------------------------------//
// VECTOR OPERATIONS
//---------------------------------------------------------------------------//

template <class T>
inline T Vector<T>::dot(const Vector<T>& x)
{
  Require(d_size == x.size());
  T val = 0;
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecDot(d_petsc_vector, const_cast<Vector<T>* >(&x)->petsc_vector(), &val);
#else
  for (int i = 0; i < d_size; i++)
    val += d_value[i] * x[i];
#endif
  return val;
}

template <class T>
inline T Vector<T>::norm(const int type)
{
  T val = 0.0;
#ifdef CALLOW_ENABLE_PETSC_OPS
  if (type == L1)
    VecNorm(d_petsc_vector, NORM_1, &val);
  else if (type == L2 or type == L2GRID)
    VecNorm(d_petsc_vector, NORM_2, &val);
  else
    VecNorm(d_petsc_vector, NORM_INFINITY, &val);
#else
  if (type == L1 or type == L1GRID)
  {
    for (int i = 0; i < d_size; i++)
      val += std::abs(d_value[i]);
  }
  else if (type == L2 or type == L2GRID)
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
  if (type == L1GRID) val /= (T)d_size;
  if (type == L2GRID) val /= std::sqrt((T)d_size);
  return val;
}

template <class T>
inline T Vector<T>::norm_residual(const Vector<T>& x, const int type)
{
  Require(d_size == x.size());
  T val = 0.0;
#ifdef CALLOW_ENABLE_PETSC_OPS
  // create the vector tmp = me - x;
  Vector<T> tmp(*this);
  VecAXPY(tmp.petsc_vector(), -1.0, const_cast<Vector<T>* >(&x)->petsc_vector());
  // and take its norm
  val = tmp.norm(type);
#else
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
#endif
  return val;
}

template <class T>
inline void Vector<T>::set(const T v)
{
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecSet(d_petsc_vector, v);
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] = v;
#endif
}

template <class T>
inline void Vector<T>::scale(const T v)
{
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecScale(d_petsc_vector, v);
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] *= v;
#endif
}

template <class T>
inline void Vector<T>::add(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecAXPY(d_petsc_vector, 1.0, const_cast<Vector<T>* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] += x[i];
#endif
}

template <class T>
inline void Vector<T>::add_a_times_x(const T a, const Vector<T>& x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecAXPY(d_petsc_vector, a, const_cast<Vector<T>* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] += a*x[i];
#endif
}

template <class T>
inline void Vector<T>::subtract(const Vector<T> &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecAXPY(d_petsc_vector, -1.0, const_cast<Vector<T>* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] -= x[i];
#endif
}

template <class T>
inline void Vector<T>::multiply(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  Vector<T> tmp(*this);
  VecPointwiseMult(d_petsc_vector, tmp.petsc_vector(), const_cast<Vector<T>* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] *= x[i];
#endif
}

template <class T>
inline void Vector<T>::divide(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  Vector<T> tmp(*this);
  VecPointwiseDivide(d_petsc_vector, tmp.petsc_vector(), const_cast<Vector<T>* >(&x)->petsc_vector());
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] /= x[i];
#endif
}

template <class T>
inline void Vector<T>::copy(const Vector &x)
{
  Require(x.size() == d_size);
#ifdef CALLOW_ENABLE_PETSC_OPS
  VecCopy(const_cast<Vector<T>* >(&x)->petsc_vector(), d_petsc_vector);
#else
  for (int i = 0; i < d_size; i++)
    d_value[i] = x[i];
#endif
}

//---------------------------------------------------------------------------//
// IO
//---------------------------------------------------------------------------//

template <class T>
void Vector<T>::display() const
{
  printf(" Vector \n");
  printf(" ---------------------------\n");
  printf("      number rows = %5i \n\n", d_size);
  if (d_size > 100)
  {
    printf("  *** vector not printed for size > 20 *** ");
    return;
  }
  for (int i = 0; i < d_size; i++)
  {
    printf(" row  %3i | %13.6e \n", i, d_value[i]);
  }
  printf("\n");
}

template <class T>
inline void Vector<T>::print_matlab(std::string filename) const
{
  FILE * f;
  f = fopen (filename.c_str(), "w");
  for (int i = 0; i < d_size; i++)
  {
    fprintf(f, "%23.16e \n", d_value[i]);
  }
  fprintf(f, "\n");
  fclose (f);
}

} // end namespace callow

#endif /* callow_VECTOR_I_HH_ */
