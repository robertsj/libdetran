//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SphericalHarmonics.cc
 *  @author Jeremy Roberts
 *  @date   Jun 30, 2011
 *  @brief  SphericalHarmonics member definitions
 */
//---------------------------------------------------------------------------//

#include "SphericalHarmonics.hh"
#include "utilities/Constants.hh"
#include "utilities/DBC.hh"
#include "utilities/GenException.hh"
#include "utilities/SoftEquivalence.hh"
#include <cmath>
#ifdef DETRAN_ENABLE_BOOST
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>
#endif

namespace detran_angle
{

// Calculate the spherical harmonic of degree l, order m
// given  cos(polar) and azimuthal angle
double SphericalHarmonics::Y_lm(const int    l,
                                const int    m,
                                const double xi,
                                const double varphi )
{
    Require (xi >= -1.0);
    Require (xi <= 1.0);
    Require (varphi >= 0.0);
    Require (varphi <= detran_utilities::two_pi);

    double sintheta = std::sqrt(1.0 - xi*xi);
    double mu  = sintheta*std::cos(varphi);
    double eta = sintheta*std::sin(varphi);
    return get_Y_lm(l, m, mu, eta, xi);

}

// Calculate the spherical harmonic of degree l, order m
// given the triplet of directional cosines
double SphericalHarmonics::Y_lm(const int    l,
                                const int    m,
                                const double mu,
                                const double eta,
                                const double xi )
{
    return get_Y_lm(l, m, mu, eta, xi);
}

// Calculate the Legendre polynomial of degree l
// given the polar angle cosine
double SphericalHarmonics::Y_lm(const int    l,
                                const double xi )
{
    return get_Y_lm(l, 0, 0, 0, xi);
}

double SphericalHarmonics::get_Y_lm(const int    l,
                                    const int    m,
                                    const double mu,
                                    const double eta,
                                    const double xi )
{
  using std::fabs;
  Require(l >= 0);
  Require(m <= l);
  Require(m >= -l);
  Require(mu >= -1.0);
  Require(mu <= 1.0);
  Require(eta >= -1.0);
  Require(eta <= 1.0);
  Require(xi >= -1.0);
  Require(xi <= 1.0);
  double unity = mu * mu + eta * eta + xi * xi;
  if (!detran_utilities::soft_equiv(1.0, unity, 1.0e-5))
  {
    Require( detran_utilities::soft_equiv( fabs(mu), 0.0 ) &&
             detran_utilities::soft_equiv( fabs(eta), 0.0 ) );
  }

  if ( l == 0 )
    return 1.0;
  else if ( l == 1 )
  {
    if ( m == -1 )
      return eta;
    else if ( m == 0 )
      return xi;
    else if ( m == 1 )
      return mu;
  }
  else if ( l == 2 )
  {
    if ( m == -2 )
      return 1.732050807568877*mu*eta;
    else if ( m == -1 )
      return 1.732050807568877*xi*eta;
    else if ( m == 0 )
      return 1.5*xi*xi-0.5;
    else if ( m == 1 )
      return 1.732050807568877*xi*mu;
    else if ( m == 2 )
      return 0.8660254037844385*(mu*mu-eta*eta);
  }
  else if ( l == 3 )
  {
    if ( m == -3 )
      return 0.7905694150420948*eta*(3*mu*mu-eta*eta);
    else if ( m == -2 )
      return 3.872983346207417*xi*mu*eta;
    else if ( m == -1 )
      return 0.6123724356957945*eta*(5.0*xi*xi-1.0);
    else if ( m == 0 )
      return 2.5*xi*xi*xi - 1.5*xi;
    else if ( m == 1 )
      return 0.6123724356957945*mu*(5.0*xi*xi-1.0);
    else if ( m == 2 )
      return 1.936491673103708*xi*(mu*mu-eta*eta);
    else if ( m == 3 )
      return 0.7905694150420948*mu*(mu*mu-3.0*eta*eta);
  }
  else
  {
#ifdef DETRAN_ENABLE_BOOST
    double phi   = std::acos(mu / std::sqrt(1.0 - xi*xi));
    double theta = std::acos(xi);
    // Normalization so that phi_00 = phi, phi_1m = current components, etc.
    // This leads to a match with definitions given by Hebert.
    double del = (m == 0) ? 1.0 : 2.0;
    double norm = std::pow(-1.0, m) *
                  std::sqrt(del * detran_utilities::four_pi / (2.*l + 1.));
    if (m >= 0)
      return norm * boost::math::spherical_harmonic_r(l, m, theta, phi);
    else
      return norm * boost::math::spherical_harmonic_i(l,-m, theta, phi);
#else
    // Degree not implemented.
    THROW("Maximum Legendre order is 3. For higher orders, enable Boost.");
	return 0;
#endif
  }
  return 0;
}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of SphericalHarmonics.cc
//---------------------------------------------------------------------------//
