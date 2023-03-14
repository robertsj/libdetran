//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Testing.cc
 *  @brief Test of TestDriver and related macros
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>

// LIST OF TEST FUNCTIONS
#define TEST_LIST                \
        FUNC(test_Testing_pass)  \
        FUNC(test_Testing_fail)

#include "TestDriver.hh"
#include <iostream>
#include <cstdlib>

using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

TEST(Testing, TestPass)
{
  EXPECT_TRUE(1==1)
}

TEST(Testing, TestFailure)
{
  EXPECT_FALSE(1==0)
}

//---------------------------------------------------------------------------//
//              end of testTesting.cc
//---------------------------------------------------------------------------//
