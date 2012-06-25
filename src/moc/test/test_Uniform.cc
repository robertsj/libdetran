//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Uniform.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of TabuchiYamamoto class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Uniform)

// Detran headers
#include "TestDriver.hh"
#include "Uniform.hh"

// Setup
/* ... */

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_Uniform()
{
  // Construct quadrature.
  Uniform q(2, 5, 13, 3, "TY");
  TEST(q.number_angles() == 4 * 5 * 3);
  TEST(q.number_angles_octant() == 5 * 3);
  TEST(q.number_azimuths_octant() == 5);
  TEST(q.number_polar_octant() == 3);
  TEST(q.number_tracks(0) == 13);
  TEST(q.number_enter(0, 0) == 2);
  TEST(q.number_enter(0, 1) == 11);
  double y = 1.0 - 0.5 / 11.0;
  TEST(soft_equiv(q.enter(0, 0).x(), 0.0));
  TEST(soft_equiv(q.enter(0, 0).y(), y));

  TEST(q.polar(13) == 1);
  TEST(q.azimuth(13) == 4);
  q.display_tracks();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
