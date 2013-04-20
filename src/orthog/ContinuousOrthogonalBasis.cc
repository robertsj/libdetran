//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ContinuousOrthogonalBasis.cc
 *  @brief  ContinuousOrthogonalBasis member definitions
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

//---------------------------------------------------------------------------//
ContinuousOrthogonalBasis::ContinuousOrthogonalBasis(const size_t order,
                                                     const vec_dbl &x,
                                                     const vec_dbl &qw)
  : OrthogonalBasis(order, x.size())
  , d_x(x)
  , d_qw(qw)
{
  // Preconditions
  Require(d_x.size() > 0);
  Require(d_x.size() == d_qw.size());

  // Allocate the weight vector, since this will store at least the qw's
  d_w = Vector::Create(d_size, 0.0);
}

//---------------------------------------------------------------------------//
ContinuousOrthogonalBasis::~ContinuousOrthogonalBasis()
{
  /* ... */
}

} // end namespace detran_orthog

//---------------------------------------------------------------------------//
//              end of file ContinuousOrthogonalBasis.cc
//---------------------------------------------------------------------------//
