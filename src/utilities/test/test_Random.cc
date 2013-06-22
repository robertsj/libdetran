//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Random.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST              \
        FUNC(test_Random)

#include "TestDriver.hh"
#include "utilities/Random.hh"
#include "utilities/SoftEquivalence.hh"
#include <cstdio>

using namespace std;
using namespace detran_test;
using namespace detran_utilities;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

double ref[] =
{ 0.33480420098231534, 0.00626434885800720, 0.80570957030638723,
  0.34539318670605823, 0.52716415307804354, 0.65810185069018345,
  0.97420811180151290, 0.62774885706890515, 0.85719300561930467,
  0.18394035853152602 };

int test_Random(int argc, char *argv[])
{

  Random R;
  for (int i = 0; i < 10; ++i)
  {
    double v = R.rnd();
    printf(" %16.10f \n", v);
    TEST(soft_equiv(v, ref[i]));
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Random.cc
//----------------------------------------------------------------------------//
