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
  : Quadrature(dim, na * np * std::pow((float)2, (int)dim), name)
  , d_number_azimuth_octant(na)
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
  size_t n = 0;
  for (size_t a = 0; a < number_azimuths_octant(); ++a)
  {
    for (size_t p = 0; p < number_polar_octant(); ++p, ++n)
    {
      d_mu[n]     = d_sin_theta[p] * d_cos_phi[a];
      d_eta[n]    = d_sin_theta[p] * d_sin_phi[a];
      d_xi[n]     = d_cos_theta[p];
      d_weight[n] = scale * d_polar_weight[p] * d_azimuth_weight[a];
      weight_tot += d_weight[n];
    } // end polar loop
  } // end azimuth loop
  weight_tot *= d_number_octants;

  verify();

  Ensure(detran_utilities::soft_equiv(weight_tot, 1.0/angular_norm(d_dimension)));
}

//---------------------------------------------------------------------------//
double ProductQuadrature::polar_weight(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_polar_weight[p];
}

//---------------------------------------------------------------------------//
double ProductQuadrature::azimuth_weight(const size_t a) const
{
  Require(a < d_number_azimuth_octant);
  return d_azimuth_weight[a];
}

//---------------------------------------------------------------------------//
void ProductQuadrature::verify() const
{
  // verify polar ordering
  double last_p = 0.0;
  for (size_t p = 0; p < d_number_polar_octant; ++p)
  {
    Insist(d_cos_theta[p] > last_p, "Non-monotonic increasing polar cosine");
    last_p = d_cos_theta[p];
    // verify azimuth
    double last_mu = 1.0;
    for (size_t a = 0; a < d_number_azimuth_octant; ++a)
    {
      Insist(d_mu[angle(a, p)] < last_mu, "Non-monotonic increasing mu");
      last_mu = d_mu[angle(a, p)];
    }
  }
}


} // end namespace detran_angle

