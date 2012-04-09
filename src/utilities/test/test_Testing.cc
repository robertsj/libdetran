//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Testing.cc
 * \author Jeremy Roberts
 * \date   Jul 14, 2011
 * \brief  Test of the simple Testing utilities.
 * \note   Copyright (C) 2011 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                \
        FUNC(test_Testing_pass)  \
        FUNC(test_Testing_fail)

// Detran headers
#include "TestDriver.hh"

// System headers
#include <iostream>
#include <cstdlib>

using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

// Test definitions.

int test_Testing_pass()
{
 TEST(1==1);
 return 0;
}

int test_Testing_fail()
{
 TESTFALSE(1==0);
 return 0;
}

//---------------------------------------------------------------------------//
//              end of testTesting.cc
//---------------------------------------------------------------------------//
