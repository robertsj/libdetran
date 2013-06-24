//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.cc
 *  @brief ProductQuadrature
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "angle/ProductQuadrature.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
ProductQuadrature::ProductQuadrature(const size_t       dim,
                                     const size_t       na,
                                     const size_t       np,
                                     const std::string &name,
                                     const bool         normalize)
  : Quadrature(dim,
               na * np * std::pow((float)2, (int)dim),
               name)
  , d_number_azimuth_octant(na)
  , d_number_polar_octant(np)
  , d_phi(d_number_azimuth_octant, 0.0)
  , d_cos_phi(d_number_azimuth_octant, 0.0)
  , d_sin_phi(d_number_azimuth_octant, 0.0)
  , d_azimuth_weight(d_number_azimuth_octant, 0.0)
  , d_cos_theta(d_number_polar_octant, 0.0)
  , d_sin_theta(d_number_polar_octant, 0.0)
  , d_polar_weight(d_number_polar_octant, 0.0)
{
  Insist(dim > 1, "Product quadratures make no sense in one dimension");
}

//----------------------------------------------------------------------------//
ProductQuadrature::~ProductQuadrature()
{
  /* ... */
}

//----------------------------------------------------------------------------//
double ProductQuadrature::phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_phi[a];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::cos_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_cos_phi[a];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::sin_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_sin_phi[a];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::azimuth_weight(const size_t a) const
{
  Require(a < d_number_azimuth_octant);
  return d_azimuth_weight[a];
}

//----------------------------------------------------------------------------//
double ProductQuadrature::cos_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_cos_theta[p];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::sin_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_sin_theta[p];
}
//----------------------------------------------------------------------------//
double ProductQuadrature::polar_weight(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_polar_weight[p];
}

//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::angle(const size_t a, const size_t p) const
{
  Require(a < d_number_azimuth_octant);
  Require(p < d_number_polar_octant);
  return p + a * d_number_polar_octant;
}
//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::azimuth(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = (angle % d_number_angles_octant) / d_number_polar_octant;
  Ensure(tmp < d_number_azimuth_octant);
  return tmp;
}
//----------------------------------------------------------------------------//
ProductQuadrature::size_t
ProductQuadrature::polar(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = angle % d_number_polar_octant;
  Ensure(tmp < d_number_polar_octant);
  return tmp;
}

//----------------------------------------------------------------------------//
void ProductQuadrature::build()
{
  // Construct the product set
  double scale = 1.0;
  if (d_dimension == 2) scale = 2.0;
  double wt_tot = 0.0;
  size_t n = 0;
  for (size_t a = 0; a < number_azimuths_octant(); ++a)
  {
    for (size_t p = 0; p < number_polar_octant(); ++p, ++n)
    {
      d_mu[n]     = d_sin_theta[p] * d_cos_phi[a];
      d_eta[n]    = d_sin_theta[p] * d_sin_phi[a];
      d_xi[n]     = d_cos_theta[p];
      d_weight[n] = scale * d_polar_weight[p] * d_azimuth_weight[a];
      wt_tot += d_weight[n];
    } // end polar loop
  } // end azimuth loop
  wt_tot *= d_number_octants;

  verify();
  Ensure(detran_utilities::soft_equiv(wt_tot, 1.0/angular_norm(d_dimension)));
}

//----------------------------------------------------------------------------//
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


} // namespace detran_angle

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.cc
//----------------------------------------------------------------------------//
