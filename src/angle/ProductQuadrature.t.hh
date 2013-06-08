//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.t.hh
 *  @brief ProductQuadrature member definitions
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_PRODUCTQUADRATURE_T_HH_
#define detran_angle_PRODUCTQUADRATURE_T_HH_

#include "ProductQuadrature.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
template <class A, class P>
ProductQuadrature<A, P>::ProductQuadrature(const size_t dim,
                                           const size_t na,
                                           const size_t np,
                                           const bool   normalize)
  : Quadrature(dim,
               na * np * std::pow((float)2, (int)dim),
               "product-" + A::name() + "-" + P::name())
  , d_number_azimuth_octant(na)
  , d_number_polar_octant(np)
{
  Insist(dim > 1, "Product quadratures make no sense in one dimension");

  // build azimuth
  A Q_A(d_number_azimuth_octant, normalize);
  d_phi = Q_A.phi();
  d_cos_phi = Q_A.cos_phi();
  d_sin_phi = Q_A.sin_phi();
  d_azimuth_weight = Q_A.weights();

  // build polar
  P Q_P(d_number_polar_octant, normalize);
  d_cos_theta = Q_P.cos_theta();
  d_sin_theta = Q_P.sin_theta();
  d_polar_weight = Q_P.weights();

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
template <class A, class P>
double ProductQuadrature<A, P>::phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_phi[a];
}
//----------------------------------------------------------------------------//
template <class A, class P>
double ProductQuadrature<A, P>::cos_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_cos_phi[a];
}
//----------------------------------------------------------------------------//
template <class A, class P>
double ProductQuadrature<A, P>::sin_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_sin_phi[a];
}
//----------------------------------------------------------------------------//
template <class A, class P>
double ProductQuadrature<A, P>::azimuth_weight(const size_t a) const
{
  Require(a < d_number_azimuth_octant);
  return d_azimuth_weight[a];
}

//----------------------------------------------------------------------------//
template <class A, class P>
double ProductQuadrature<A, P>::cos_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_cos_theta[p];
}
//----------------------------------------------------------------------------//
template <class A, class P>
double ProductQuadrature<A, P>::sin_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_sin_theta[p];
}
//----------------------------------------------------------------------------//
template <class A, class P>
double ProductQuadrature<A, P>::polar_weight(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_polar_weight[p];
}

//----------------------------------------------------------------------------//
template <class A, class P>
typename ProductQuadrature<A, P>::size_t
ProductQuadrature<A, P>::angle(const size_t a, const size_t p) const
{
  Require(a < d_number_azimuth_octant);
  Require(p < d_number_polar_octant);
  return p + a * d_number_polar_octant;
}
//----------------------------------------------------------------------------//
template <class A, class P>
typename ProductQuadrature<A, P>::size_t
ProductQuadrature<A, P>::azimuth(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = (angle % d_number_angles_octant) / d_number_polar_octant;
  Ensure(tmp < d_number_azimuth_octant);
  return tmp;
}
//----------------------------------------------------------------------------//
template <class A, class P>
typename ProductQuadrature<A, P>::size_t
ProductQuadrature<A, P>::polar(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = angle % d_number_polar_octant;
  Ensure(tmp < d_number_polar_octant);
  return tmp;
}

//----------------------------------------------------------------------------//
template <class A, class P>
void ProductQuadrature<A, P>::verify() const
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

#endif /* detran_angle_PRODUCTQUADRATURE_T_HH_ */

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.t.hh
//----------------------------------------------------------------------------//
