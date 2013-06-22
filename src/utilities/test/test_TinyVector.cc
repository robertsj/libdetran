//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_TinyVector.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST              \
        FUNC(test_TinyVector)

#include "TestDriver.hh"
#include "utilities/TinyVector.hh"

using namespace detran_test;
using namespace detran_utilities;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_TinyVector(int argc, char *argv[])
{
  TinyVector<int, 2> V2(1, 2);
  TEST(V2.size() == 2);
  TEST(V2[0] == 1);
  TEST(V2[1] == 2);

  TinyVector<int, 3> V3(1, 2, 3);
  TEST(V3.size() == 3);
  TEST(V3[0] == 1);
  TEST(V3[1] == 2);
  TEST(V3[2] == 3);

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_TinyVector.cc
//----------------------------------------------------------------------------//
