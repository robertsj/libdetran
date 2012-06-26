//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Collocated.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Collocated class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Collocated)

// Detran headers
#include "TestDriver.hh"
#include "Collocated.hh"

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

int test_Collocated()
{
  // Construct quadrature.
  Collocated q(2, 5, 1, 2, "TY");
  TEST(q.number_angles()          == 4 * 5 * 2);
  TEST(q.number_angles_octant()   == 5 * 2);
  TEST(q.number_azimuths_octant() == 5);
  TEST(q.number_polar_octant()    == 2);
  TEST(q.number_tracks(0)         == 10);
  TEST(q.number_enter(0, 0)       == 1);
  TEST(q.number_enter(0, 1)       == 9);
  TEST(soft_equiv(q.enter(0, 0).x(), 0.0));
  TEST(soft_equiv(q.enter(0, 0).y(), 9.444444444444444e-01));

//  TEST(q.polar(13) == 1);
//  TEST(q.azimuth(13) == 4);
//  q.display();
  //q.display_tracks();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
