//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Iterators.i.hh
 *  @brief Iterator inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_utilities_ITERATORS_I_HH_
#define detran_utilities_ITERATORS_I_HH_

namespace detran_utilities
{

//----------------------------------------------------------------------------//
template <class C>
Reversible<C>::Reversible(container_type *c, bool d)
  : d_container(c)
{
  if (d)
  {
    d_p = d_container->begin();
    d_increment = 1;
  }
  else
  {
    d_p = d_container->end() - 1;
    d_increment = -1;
  }
}

//----------------------------------------------------------------------------//
template <class C>
Reversible<C>& Reversible<C>::operator+(const int n)
{
  d_p += d_increment * n;
  return *this;
}

//----------------------------------------------------------------------------//
template <class C>
Reversible<C>& Reversible<C>::operator++()
{
  d_p += d_increment;
  return *this;
}

//----------------------------------------------------------------------------//
template <class C>
Reversible<C>& Reversible<C>::operator++(int)
{
  d_p += d_increment;
  return *this;
}

//----------------------------------------------------------------------------//
template <class C>
bool Reversible<C>::operator==(const Reversible& o) const
{
  return (d_p == o.d_p);
}

//----------------------------------------------------------------------------//
template <class C>
bool Reversible<C>::operator!=(const Reversible& o) const
{
  return !(*this == o);
}

//----------------------------------------------------------------------------//
template <class C>
typename Reversible<C>::value_type& Reversible<C>::operator*()
{
  return *d_p;
}

//----------------------------------------------------------------------------//
template <class C>
typename Reversible<C>::base_iterator& Reversible<C>::operator->()
{
  return d_p;
}

} // end namespace detran_utilities

#endif /* detran_utilities_ITERATORS_I_HH_ */

//----------------------------------------------------------------------------//
//              end of file Iterators.hh
//----------------------------------------------------------------------------//
