//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Testing.cc
 *  @brief Test of TestDriver and related macros
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>

using namespace std;

TEST(Testing, TestPass)
{
  EXPECT_TRUE(1==1);
}

TEST(Testing, TestFailure)
{
  EXPECT_FALSE(1==0);
}

//---------------------------------------------------------------------------//
//              end of testTesting.cc
//---------------------------------------------------------------------------//
