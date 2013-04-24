//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   OrthogonalBasis.cc
 *  @brief  OrthogonalBasis member definitions
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#include "orthog/OrthogonalBasis.hh"
#include <cmath>

namespace detran_orthog
{

//---------------------------------------------------------------------------//
OrthogonalBasis::OrthogonalBasis(const size_t order,
                                 const size_t size,
                                 const bool   orthonormal)
  : d_order(order)
  , d_size(size)
  , d_orthonormal(orthonormal)
{
  Require(d_order < d_size);
}

//---------------------------------------------------------------------------//
OrthogonalBasis::~OrthogonalBasis()
{
  /* ... */
}

//---------------------------------------------------------------------------//
void OrthogonalBasis::compute_a()
{
  Require(d_basis);

  if (!d_a) d_a = Vector::Create(d_order + 1);
  d_a->set(1.0);

  // If we're orthonormal, then we pre-normalize the basis.
  for (size_t i = 0; i <= d_order; ++i)
  {
    Vector row(d_size, &(*d_basis)(i, 0));
    Vector tmp(row);
    if (d_w)
    {
      tmp.multiply(d_w);
    }
    double val = 1.0 / std::sqrt(tmp.dot(row));

    if (d_orthonormal)
    {
      row.scale(val);
    }
    else
    {
      (*d_a)[i] = val;
    }
  }

}

//---------------------------------------------------------------------------//
template ORTHOG_EXPORT class detran_utilities::SP<OrthogonalBasis>;

} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file OrthogonalBasis.cc
//---------------------------------------------------------------------------//
