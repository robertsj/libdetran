//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SphericalHarmonics.cc
 * \author Jeremy Roberts
 * \date   Jun 30, 2011
 * \brief
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

// Detran
#include "SphericalHarmonics.hh"

// Utilities
#include "Constants.hh"
#include "DBC.hh"
#include "GenException.hh"
#include "SoftEquivalence.hh"

// System
#include <cmath>

namespace detran
{

// Calculate the spherical harmonic of degree l, order m
// given  cos(polar) and azimuthal angle
double SphericalHarmonics::Y_lm(int l,
                                int m,
                                double xi,
                                double varphi )
{

    using std::sqrt;

    Require (xi >= -1.0);
    Require (xi <= 1.0);
    Require (varphi >= 0.0);
    Require (varphi <= two_pi);

    double sintheta = sqrt(1.0 - xi*xi);
    double mu     = sintheta*cos(varphi);
    double eta    = sintheta*sin(varphi);
    return get_Y_lm(l,m,mu,eta,xi);

}

// Calculate the spherical harmonic of degree l, order m
// given the triplet of directional cosines
double SphericalHarmonics::Y_lm( int l,
                              int m,
                              double mu,
                              double eta,
                              double xi )
{

    return get_Y_lm(l,m,mu,eta,xi);

}

// Calculate the Legendre polynomial of degree l
// given the polar angle cosine
double SphericalHarmonics::Y_lm( int l,
                              double xi )
{

    return get_Y_lm(l,0,0,0,xi);

}

double  SphericalHarmonics::get_Y_lm(int l,
                                     int m,
                                     double mu,
                                     double eta,
                                     double xi )
{
  using std::fabs;
  Require (l >= 0);
  Require (m <= l);
  Require (m >= -l);
  Require ( (mu >= -1.0) && (mu <= 1.0) );
  Require ( (eta >= -1.0) && (eta <= 1.0) );
  Require ( (xi >= -1.0) && (xi <= 1.0) );
  double unity = mu * mu + eta * eta + xi * xi;
  if (!soft_equiv(1.0, unity, 1.0e-5))
  {
    Require( soft_equiv( fabs(mu), 0.0 ) &&
             soft_equiv( fabs(eta), 0.0 ) );
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

  // Degree not implemented.
  THROW("l too large; maximum Legendre order is 3.");

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of SphericalHarmonics.cc
//---------------------------------------------------------------------------//
