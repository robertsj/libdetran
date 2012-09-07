//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_BoundarySN.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of BoundarySN class
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_BoundarySN)

// Detran headers
#include "TestDriver.hh"
#include "boundary/BoundarySN.hh"
#include "angle/GaussLegendre.hh"
#include "angle/QuadrupleRangle.hh"

// Setup
#include "quadrature_fixture.hh"

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

int test_BoundarySN_basic(int argc, char *argv[])
{
  // Get quadrature fixture
  SP_quadrature q = gausslegendre_fixture();
  TEST(q);
  TEST(q->number_angles()  == 8);
  TEST(q->number_octants() == 2);
  TEST(q->number_angles_octant() == 4);
  TEST(soft_equiv(q->mu(0, 0),  0.9602898564975));
  TEST(soft_equiv(q->mu(1, 0), -0.9602898564975));
  TEST(soft_equiv(q->weight(0), 0.1012285362904));
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_BoundarySN.cc
//---------------------------------------------------------------------------//
