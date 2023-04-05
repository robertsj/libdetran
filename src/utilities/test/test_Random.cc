//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Random.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "utilities/Random.hh"
#include "utilities/SoftEquivalence.hh"
#include <cstdio>

using namespace std;
using namespace detran_utilities;

TEST(Random, Basic)
{
  const double expected_vals[] =
  { 0.33480420098231534, 0.00626434885800720, 0.80570957030638723,
    0.34539318670605823, 0.52716415307804354, 0.65810185069018345,
    0.97420811180151290, 0.62774885706890515, 0.85719300561930467,
    0.18394035853152602 };

  Random R;
  for (int i = 0; i < 10; ++i)
  {
    double got = R.rnd();
    double expect = expected_vals[i];
    printf(" %16.10f \n", got);
    EXPECT_NEAR(got, expect, 1.0e-12);
  }
}

//----------------------------------------------------------------------------//
//              end of test_Random.cc
//----------------------------------------------------------------------------//
