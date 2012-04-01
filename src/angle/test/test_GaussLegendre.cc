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
#define TEST_LIST                     \
        FUNC(test_GaussLegendre_basic)

// Detran headers
#include "TestDriver.hh"
#include "GaussLegendre.hh"

// Setup
#include "quadrature_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace detran_utils;
using namespace std;

int main(int argc, char *argv[])
{
  // Initialize tests.
  int test = TestDriver::initialize(argc, argv);
  // Perform the test and return result.
  return TestDriver::evaluate((*test_table[test])());
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_GaussLegendre_basic()
{
  // Get quadrature fixture
  SP_quadrature q = gausslegendre_fixture();

  // Finish me.

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_GaussLegendre.cc
//---------------------------------------------------------------------------//
