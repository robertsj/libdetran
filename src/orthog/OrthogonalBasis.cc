//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  OrthogonalBasis.cc
 *  @brief OrthogonalBasis member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "orthog/OrthogonalBasis.hh"
#include <cmath>

namespace detran_orthog
{

//----------------------------------------------------------------------------//
OrthogonalBasis::OrthogonalBasis(const Parameters &p)
  : d_order(p.order)
  , d_size(p.size)
  , d_orthonormal(p.orthonormal)
  , d_even_only(p.even_only)
{
  Requirev(d_order < d_size, AsString(d_order)+" < "+AsString(d_size));
}

//----------------------------------------------------------------------------//
OrthogonalBasis::~OrthogonalBasis()
{
  /* ... */
}

//----------------------------------------------------------------------------//
void OrthogonalBasis::compute_a()
{
  Require(d_basis);

  if (!d_a) d_a = Vector::Create(d_order + 1);
  d_a->set(1.0);

  // If we're orthonormal, then we pre-normalize the basis.
  for (size_t i = 0; i <= d_order; ++i)
  {
    Vector column(d_size, &(*d_basis)(0, i));
    Vector tmp(column);
    if (d_w)
    {
      tmp.multiply(d_w);
    }
    double val = 1.0 / tmp.dot(column);

    if (d_orthonormal)
    {
      // Then scale V so that V'V = 1
      column.scale(std::sqrt(val));
    }
    else
    {
      // Scale so that value * V'V = 1
      (*d_a)[i] = val;
    }
  }

}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file OrthogonalBasis.cc
//----------------------------------------------------------------------------//
