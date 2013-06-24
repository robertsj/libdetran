//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Iterators.hh
 *  @brief Iterators class definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 *
 *  The iterators defined are essentially convenience wrappers around
 *  STL iterators.  While Boost could be used, doing so would entail
 *  a lot of extra stuff, leading to longer compilation times that we
 *  just don't need.
 */
//----------------------------------------------------------------------------//

#ifndef detran_utilities_ITERATORS_HH_
#define detran_utilities_ITERATORS_HH_

namespace detran_utilities
{

/// Lightweight wrapper around a STL iterator to enable reverse in-place
template <class C>
class Reversible
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef typename C::value_type      value_type;
  typedef C                           container_type;
  typedef typename C::iterator        base_iterator;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Constructor with container pointer and flag indicating forward progress
  Reversible(container_type *c, bool d = true);
  /// Increment by an arbitrary number
  Reversible& operator+(const int n);
  /// Increment by one
  Reversible& operator++();
  /// Increment by one (for post-fix)
  Reversible& operator++(int);
  //@{
  /// Equivalence
  bool operator==(const Reversible& o) const;
  bool operator!=(const Reversible& o) const;
  //@}
  /// Dereference
  value_type& operator*();
  base_iterator& operator->();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Base container
  container_type* d_container;
  /// Base iterator
  base_iterator d_p;
  /// Increment (+1 or -1)
  int d_increment;

};

} // end namespace detran_utilities

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "Iterators.i.hh"

#endif /* detran_utilities_ITERATORS_HH_ */

//----------------------------------------------------------------------------//
//              end of file Iterators.hh
//----------------------------------------------------------------------------//
