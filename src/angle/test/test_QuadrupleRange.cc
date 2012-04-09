//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_QuadrupleRange.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of QuadrupleRange class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_QuadrupleRange_basic)

// Detran headers
#include "TestDriver.hh"
#include "QuadrupleRange.hh"

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

int test_QuadrupleRange_basic()
{
  // Get quadrature fixture
  SP_quadrature q = quadruplerange_fixture();
  q->display();
  TEST(q);
  TEST(q->number_angles()   == 8);
  TEST(q->number_octants()  == 4);
  TEST(soft_equiv(q->mu(0, 0),  0.2582870761957));
  TEST(soft_equiv(q->mu(0, 1),  0.7417129238043));
  TEST(soft_equiv(q->eta(0, 0), 0.7417129238043));
  TEST(soft_equiv(q->eta(0, 1), 0.2582870761957));
  TEST(soft_equiv(q->weight(0), 1.5707963267949));
  TEST(soft_equiv(q->weight(1), 1.5707963267949));
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_QuadrupleRange.cc
//---------------------------------------------------------------------------//
