//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Segment.cc
 *  @brief Test of Segment class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST          \
        FUNC(test_Segment)

#include "TestDriver.hh"
#include "geometry/Segment.hh"

using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_Segment(int argc, char *argv[])
{
  Segment s(0, 1.0);
  TEST(s.region() == 0);
  TEST(soft_equiv(s.length(), 1.0));
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Segment.cc
//----------------------------------------------------------------------------//
