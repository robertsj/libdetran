//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ProductQuadrature.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  ProductQuadrature class definition.
 */
//---------------------------------------------------------------------------//

#include "ProductQuadrature.hh"
#include <cmath>

namespace detran_angle
{

ProductQuadrature::ProductQuadrature(const size_t dim,
                                     const size_t na,
                                     const size_t np,
                                     std::string  name)
  : Quadrature(dim, na * np * std::pow(2, dim), name)
  , d_number_azimuths_octant(na)
  , d_number_polar_octant(np)
  , d_phi(na, 0.0)
  , d_cos_phi(na, 0.0)
  , d_sin_phi(na, 0.0)
  , d_azimuth_weight(na, 0.0)
  , d_cos_theta(np, 0.0)
  , d_sin_theta(np, 0.0)
  , d_polar_weight(np, 0.0)
{
  // Preconditions
  Insist(dim > 1, "Product quadratures make no sense in one dimension");
}

ProductQuadrature::~ProductQuadrature()
{
  /* ... */
}

} // end namespace detran_angle



