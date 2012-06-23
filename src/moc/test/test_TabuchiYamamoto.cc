//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_TabuchiYamamoto.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of TabuchiYamamoto class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_TabuchiYamamoto)

// Detran headers
#include "TestDriver.hh"
#include "TabuchiYamamoto.hh"

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

int test_TabuchiYamamoto()
{
  // Construct quadrature.
  TabuchiYamamoto q(1);
  TEST(q.number_polar() == 1);
  TEST(soft_equiv(q.sin_theta(0),  0.798184));
  TEST(soft_equiv(q.weight(0),     1.000000));
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_TabuchiYamamoto.cc
//---------------------------------------------------------------------------//
