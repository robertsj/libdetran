//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ProductQuadrature.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  ProductQuadrature class definition.
 */
//---------------------------------------------------------------------------//

#include "ProductQuadrature.hh"
#include "utilities/SoftEquivalence.hh"
#include <cmath>

namespace detran_angle
{

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
ProductQuadrature::~ProductQuadrature()
{
  /* ... */
}

//---------------------------------------------------------------------------//
void ProductQuadrature::build_product_quadrature()
{
  double scale = 1.0;
  if (d_dimension == 2) scale = 2.0;
  double weight_tot = 0.0;
  int n = 0;
  for (int a = 0; a < number_azimuths_octant(); ++a)
  {
    for (int p = 0; p < number_polar_octant(); ++p, ++n)
    {
      d_mu[n]     = d_sin_theta[p] * d_cos_phi[a];
      d_eta[n]    = d_sin_theta[p] * d_sin_phi[a];
      d_xi[n]     = d_cos_theta[p];
      d_weight[n] = scale * d_polar_weight[p] * d_azimuth_weight[a];
      weight_tot += d_weight[n];
    } // end polar loop
  } // end azimuth loop
  weight_tot *= d_number_octants;

  // Postconditions
  Ensure(detran_utilities::soft_equiv(weight_tot, 1.0/angular_norm(d_dimension)));
}

} // end namespace detran_angle

