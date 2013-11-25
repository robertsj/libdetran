//----------------------------------*-C++-*---------------------------0-------//
/**
 *  @file  TinyVector.i.hh
 *  @brief TinyVector inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//
#ifndef detran_utilities_TINYVECTOR_I_HH_
#define detran_utilities_TINYVECTOR_I_HH_

#include "utilities/DBC.hh"
#include "utilities/TinyVector.hh"
#include <cstring>

namespace detran_utilities
{

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector()
{
  for (size_t i = 0; i < N; ++i)
    d_values[i] = T();
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(const TinyVector &V)
{
  std::memcpy(&d_values[0], &V[0], N * sizeof(value_type));
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(const value_type *values)
{
  std::memcpy(&d_values[0], &values[0], N * sizeof(value_type));
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(cvt v0)
{
  StaticAssert(N == 1);
  d_values[0] = v0;
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(cvt v0, cvt v1)
{
  StaticAssert(N == 2);
  d_values[0] = v0;
  d_values[1] = v1;
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(cvt v0, cvt v1, cvt v2)
{
  StaticAssert(N == 3);
  d_values[0] = v0;
  d_values[1] = v1;
  d_values[2] = v2;
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(cvt v0, cvt v1, cvt v2, cvt v3)
{
  StaticAssert(N == 4);
  d_values[0] = v0;
  d_values[1] = v1;
  d_values[2] = v2;
  d_values[3] = v3;
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(cvt v0, cvt v1, cvt v2, cvt v3, cvt v4)
{
  StaticAssert(N == 5);
  d_values[0] = v0;
  d_values[1] = v1;
  d_values[2] = v2;
  d_values[3] = v3;
  d_values[4] = v4;
}

//----------------------------------------------------------------------------//
template<class T, int N>
TinyVector<T, N>::TinyVector(cvt v0, cvt v1, cvt v2, cvt v3, cvt v4, cvt v5)
{
  StaticAssert(N == 6);
  d_values[0] = v0;
  d_values[1] = v1;
  d_values[2] = v2;
  d_values[3] = v3;
  d_values[4] = v4;
  d_values[5] = v5;
}


//----------------------------------------------------------------------------//
template<class T, int N>
const T& TinyVector<T, N>::operator[](const size_type n) const
{
  Require(n < N);
  return d_values[n];
}

//----------------------------------------------------------------------------//
template<class T, int N>
T& TinyVector<T, N>::operator[](const size_type n)
{
  Require(n < N);
  return d_values[n];
}

//----------------------------------------------------------------------------//
template<class T, int N>
typename TinyVector<T, N>::iterator TinyVector<T, N>::begin()
{
  return &d_values[0];
}

//----------------------------------------------------------------------------//
template<class T, int N>
typename TinyVector<T, N>::iterator TinyVector<T, N>::end()
{
  return &d_values[N];
}

//----------------------------------------------------------------------------//
template<class T, int N>
typename TinyVector<T, N>::const_iterator TinyVector<T, N>::begin() const
{
  return &d_values[0];
}

//----------------------------------------------------------------------------//
template<class T, int N>
typename TinyVector<T, N>::const_iterator TinyVector<T, N>::end() const
{
  return &d_values[N];
}

//----------------------------------------------------------------------------//
template<class T, int N>
typename TinyVector<T, N>::size_type TinyVector<T, N>::size() const
{
  return N;
}

} // end namespace detran_utilities

#endif /* detran_utilities_TINYVECTOR_I_HH_ */

//----------------------------------------------------------------------------//
//              end of file TinyVector.i.hh
//----------------------------------------------------------------------------//
