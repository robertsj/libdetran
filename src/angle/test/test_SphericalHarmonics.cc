//----------------------------------*-C++-*-------------------=---------------//
/**
 *  @file  test_SphericalHarmonics.cc
 *  @brief Test of SphericalHarmonics class.
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                                 \
        FUNC(test_SphericalHarmonics)             \
        FUNC(test_SphericalHarmonics_integration)

#include "utilities/TestDriver.hh"
#include "angle/SphericalHarmonics.hh"
#include "angle/QuadratureFactory.hh"
#ifdef DETRAN_ENABLE_BOOST
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>
using boost::math::spherical_harmonic_r;
using boost::math::spherical_harmonic_i;
using boost::math::factorial;
#endif
#include "utilities/Constants.hh"
#include <cstdio>
// Setup
/* ... */

using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using std::cout;
using std::endl;
using std::printf;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_SphericalHarmonics(int argc, char *argv[])
{
  double mu  = 0.350021174581540677777041;
  double eta = 0.350021174581540677777041;
  double xi  = 0.868890300722201205229788;
  TEST(soft_equiv(SphericalHarmonics::Y_lm(0, 0, mu, eta, xi), 1.0));
  TEST(soft_equiv(SphericalHarmonics::Y_lm(1,-1, mu, eta, xi), eta));
  TEST(soft_equiv(SphericalHarmonics::Y_lm(1, 0, mu, eta, xi), xi));
  TEST(soft_equiv(SphericalHarmonics::Y_lm(1, 1, mu, eta, xi), mu));
  return 0;
}

/*
 *  This test compares the accuracy of 3-D quadratures.  In
 *  particular, we test a quadrature on the following
 *  moments:
 *
 *    (2/pi)int( cos(phi)^m * sin(phi)^(l+1) * sin(theta)^l, phi=0..pi/2, theta=0..pi/2)
 *
 *    l, m   even
 *
 *  = 1/(m+1)                                 for l = 0
 *    (l/2-1/2)!(m/2-1/2)! / (l/2+m/2+1/2)!   for l >= 0
 *
 *  The level symmetric quadrature integrates these exactly for
 *  N >= k + l, where N is the order (= twice the polar/octant)
 *
 *  This could be made a public function at some point.
 */
int test_SphericalHarmonics_integration(int argc, char *argv[])
{
  InputDB::SP_input inp(new InputDB());
  inp->put<int>("quad_number_polar_octant",   4);
  inp->put<int>("quad_number_azimuth_octant", 4);
  inp->put<std::string>("quad_type", "asqr-asdr");
  QuadratureFactory::SP_quadrature Q = QuadratureFactory::build(inp, 3);

  int L = 6;

  vec2_dbl err(L + 1, vec_dbl(L + 1, 0.0));

  for (int ll = 0; ll <= L; ++ll)
  {
    int l = 2*ll;
    for (int mm = 0; mm <= L; ++mm)
    {
      int m = 2*mm;

      // Compute the reference
      double ref;
      if (l == 0)
      {
        ref = 1.0 / (m + 1.0);
      }
      else
      {
        ref = 1.0 / (l + m + 1.0);
        for (int i = 1; i < l; i+=2)
        {
          double den = m + i;
          ref *= (double)i / den;
        }
      }
      ref *= four_pi;

      double val = 0.0;
      for (int o = 0; o < 8; ++o)
      {
        for (int a = 0; a < Q->number_angles_octant(); ++a)
        {
          double mu  = Q->mu(o, a);
          //double eta = Q->eta(o, a);
          double xi  = Q->xi(o, a);
          val += Q->weight(a) * std::pow(mu, m) * std::pow(xi, l);
        }
      }
      err[ll][mm] = 100.0*(val-ref)/ref;
    }
  }

  // Print a table of errors
  cout << " NUMBER OF ANGLES = " << Q->number_angles() << endl;
  printf("  l \\ m ");
  for (int mm = 0; mm <= L; ++mm)
    printf("  %3i         ", 2*mm);
  printf("\n ------");
  for (int mm = 0; mm <= L; ++mm)
    printf("--------------");
  printf("\n");
  for (int ll = 0; ll <= L; ++ll)
  {
    int l = 2*ll;
    printf(" %3i |", l);
    for (int mm = 0; mm <= L; ++mm)
    {
      int m = 2*mm;
      printf(" %12.6f ", err[ll][mm]);
    }
    printf("\n");
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_SphericalHarmonics.cc
//----------------------------------------------------------------------------//
