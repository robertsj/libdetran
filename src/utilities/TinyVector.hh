//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TinyVector.hh
 *  @brief TinyVector class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_utilities_TINYVECTOR_HH_
#define detran_utilities_TINYVECTOR_HH_

#include <iterator>

namespace detran_utilities
{

/**
 *  @class TinyVector
 *  @brief Small, compile-time N-vector
 */
template <class T, int N>
class TinyVector
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef T                                       value_type;
  typedef unsigned int                            size_type;
  typedef value_type&                             reference;
  typedef const value_type&                       const_reference;
  typedef value_type*                             pointer;
  typedef const value_type*                       const_pointer;
  typedef pointer                                 iterator;
  typedef const_pointer                           const_iterator;
  typedef std::reverse_iterator<iterator>         reverse_iterator;
  typedef std::reverse_iterator<const_iterator>   const_reverse_iterator;
  typedef const value_type                        cvt;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Default constructor
  TinyVector();
  /// Copy constructor
  TinyVector(const TinyVector &V);
  /// Constructor from values
  TinyVector(const value_type *values);
  //@{
  /// Construction from explicit values
  TinyVector(cvt v0);
  TinyVector(cvt v0, cvt v1);
  TinyVector(cvt v0, cvt v1, cvt v2);
  TinyVector(cvt v0, cvt v1, cvt v2, cvt v3);
  TinyVector(cvt v0, cvt v1, cvt v2, cvt v3, cvt v4);
  TinyVector(cvt v0, cvt v1, cvt v2, cvt v3, cvt v4, cvt v5);
  TinyVector(cvt v0, cvt v1, cvt v2, cvt v3, cvt v4, cvt v5, cvt v6);
  //@}

  //@{
  /// Element access
  const T& operator[](const size_type n) const;
  T& operator[](const size_type n);
  //@}

  //@{
  /// Iterator access
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  //@}

  /// Vector size
  size_type size() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Values
  T d_values[N];

};

} // end namespace detran_utilities

//----------------------------------------------------------------------------//
// INLINE MEMBERS
//----------------------------------------------------------------------------//

#include "TinyVector.i.hh"

#endif /* detran_utilities_TINYVECTOR_HH_ */

//----------------------------------------------------------------------------//
//              end of file TinyVector.hh
//----------------------------------------------------------------------------//
