//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_SphericalHarmonics.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of SphericalHarmonics class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_SphericalHarmonics)

// Detran headers
#include "TestDriver.hh"
#include "SphericalHarmonics.hh"

// Setup
/* ... */

using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_SphericalHarmonics(int argc, char *argv[])
{
  double val;
  double mu  = 0.350021174581540677777041;
  double eta = 0.350021174581540677777041;
  double xi  = 0.868890300722201205229788;
  TEST(soft_equiv(SphericalHarmonics::Y_lm(0, 0, mu, eta, xi), 1.0));
  TEST(soft_equiv(SphericalHarmonics::Y_lm(1,-1, mu, eta, xi), eta));
  TEST(soft_equiv(SphericalHarmonics::Y_lm(1, 0, mu, eta, xi), xi));
  TEST(soft_equiv(SphericalHarmonics::Y_lm(1, 1, mu, eta, xi), mu));
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_SphericalHarmonics.cc
//---------------------------------------------------------------------------//
