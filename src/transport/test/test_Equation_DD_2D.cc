//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Equation_DD_2D.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of Equation_DD_2D.cc
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//


// LIST OF TEST FUNCTIONS
#define TEST_LIST                  \
        FUNC(test_Equation_DD_2D)

// Detran headers
#include "TestDriver.hh"

// Setup
#include "input_fixture.hh"
#include "mesh_fixture.hh"
#include "material_fixture.hh"

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

// Test definitions.
int test_Equation_DD_2D()
{



}

//---------------------------------------------------------------------------//
//              end of test_Equation_DD_2D.cc
//---------------------------------------------------------------------------//
