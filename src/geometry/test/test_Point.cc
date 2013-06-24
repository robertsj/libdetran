//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_Point.cc
 *  @author Jeremy Roberts
 *  @note   Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_Point)

#include "TestDriver.hh"
#include "Point.hh"

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

int test_Point(int argc, char *argv[])
{
  Point p1(1.0, 2.0);
  TEST(soft_equiv(p1.x(), 1.0));
  TEST(soft_equiv(p1.y(), 2.0));
  Point p2 =  2.0 * p1;
  TEST(soft_equiv(p2.x(), 2.0));
  TEST(soft_equiv(p2.y(), 4.0));
  Point p3 = p1 + p2;
  TEST(soft_equiv(p3.x(), 3.0));
  TEST(soft_equiv(p3.y(), 6.0));
  Point p4 = p1 - p2;
  TEST(soft_equiv(p4.x(), -1.0));
  TEST(soft_equiv(p4.y(), -2.0));
  TEST(soft_equiv(distance(p1, p2), 2.23606797749979));
  //
  p1 = Point(0.75, 0.00);
  p2 = Point(0.50, 0.50);
  TEST(soft_equiv(distance(p2, p1), 0.559016994374947));
  TEST(soft_equiv(distance(p1, p2), 0.559016994374947));
  p1 = p2 + 0.5;
  TEST(soft_equiv(p1.x(), 1.0));
  TEST(soft_equiv(p1.y(), 1.0));
  TEST(soft_equiv(p1.z(), 0.5));
  p2 = p1 - 0.1;
  TEST(soft_equiv(p2.x(), 0.9));
  TEST(soft_equiv(p2.y(), 0.9));
  TEST(soft_equiv(p2.z(), 0.4));

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Point.cc
//---------------------------------------------------------------------------//
