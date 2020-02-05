//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ContinuousOrthogonalBasis.cc
 *  @brief ContinuousOrthogonalBasis member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

//----------------------------------------------------------------------------//
ContinuousOrthogonalBasis::ContinuousOrthogonalBasis(const Parameters &p)
  : OrthogonalBasis(p)
  , d_x(p.x)
  , d_qw(p.qw)
  , d_lower_bound(p.lower_bound)
  , d_upper_bound(p.upper_bound)
{
  Require(d_x.size() > 0);
  Require(d_x.size() == d_qw.size());
  Require(d_lower_bound < d_upper_bound);
  Require(d_lower_bound <= d_x[0]);
  Require(d_upper_bound >= d_x[d_x.size()-1]);

  // Allocate the weight vector, since this will store at least the qw's
  d_w = Vector::Create(d_size, 0.0);
}

//----------------------------------------------------------------------------//
ContinuousOrthogonalBasis::~ContinuousOrthogonalBasis()
{
  /* ... */
}

} // end namespace detran_orthog

//----------------------------------------------------------------------------//
//              end of file ContinuousOrthogonalBasis.cc
//----------------------------------------------------------------------------//
