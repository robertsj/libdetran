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
  // Uniform q(2, 9, 20, 3, "TY");
  Uniform q(2, 1, 3, 1, "TY");

//  TEST(q.number_angles() == 4 * 9 * 3);
//  TEST(q.number_angles_octant() == 9 * 3);
//  TEST(q.number_azimuths_octant() == 9);
//  TEST(q.number_polar() == 3);
//  TEST(q.number_tracks(0) == 20);
//  TEST(q.number_enter(0, 0) == 2);
//  TEST(q.number_enter(0, 1) == 18);
//  double y = 1.0 - 0.5 / 18.0;
//  TEST(soft_equiv(q.enter(0, 0).x(), 0.0));
//  TEST(soft_equiv(q.enter(0, 0).y(), y));
  q.display();
  q.display_tracks();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
