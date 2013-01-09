//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   OrthogonalBasis.cc
 *  @brief  OrthogonalBasis member definitions
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

//---------------------------------------------------------------------------//
OrthogonalBasis::OrthogonalBasis(const size_t order, const size_t size)
  : d_order(order)
  , d_size(size)
{
  // Preconditions
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
  // Preconditions
  Require(d_basis);

  if (!d_a) d_a = Vector::Create(d_order + 1, 0.0);

  for (size_t i = 0; i <= d_order; ++i)
  {
    Vector row(d_size, &(*d_basis)(i, 0));
    double den = row.norm(callow::L2);
    den *= den;
    (*d_a)[i] = 1.0 / den;
  }

}

} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file OrthogonalBasis.cc
//---------------------------------------------------------------------------//
