//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_GaussLegendre.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of GaussLegendre class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_GaussLegendre_basic)

// Detran headers
#include "TestDriver.hh"
#include "GaussLegendre.hh"

// Setup
#include "quadrature_fixture.hh"

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

int test_GaussLegendre_basic(int argc, char *argv[])
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
//              end of test_GaussLegendre.cc
//---------------------------------------------------------------------------//
