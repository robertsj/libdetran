//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   test_Point.cc
 *  @author Jeremy Roberts
 *  @note   Copyright (C) 2012 Jeremy Roberts.
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "Point.hh"

using namespace detran_geometry;

TEST(Point, PointAccessors)
{
  Point p1(1.0, 2.0);
  EXPECT_EQ(p1.x(), 1.0);
  EXPECT_EQ(p1.y(), 2.0);
  EXPECT_EQ(p1.z(), 0.0);
}

TEST(Point, PointOperators)
{
  Point p1(1.0, 2.0);

  Point p2 =  2.0 * p1;
  EXPECT_DOUBLE_EQ(p2.x(), 2.0);
  EXPECT_DOUBLE_EQ(p2.y(), 4.0);
  EXPECT_DOUBLE_EQ(p2.x(), 2.0);

  Point p3 = p1 + p2;
  EXPECT_DOUBLE_EQ(p3.x(), 3.0);
  EXPECT_DOUBLE_EQ(p3.y(), 6.0);
  Point p4 = p1 - p2;
  EXPECT_DOUBLE_EQ(p4.x(), -1.0);
  EXPECT_DOUBLE_EQ(p4.y(), -2.0);
  EXPECT_DOUBLE_EQ(distance(p1, p2), 2.23606797749979);
  //
  p1 = Point(0.75, 0.00);
  p2 = Point(0.50, 0.50);
  EXPECT_DOUBLE_EQ(distance(p2, p1), 0.559016994374947);
  EXPECT_DOUBLE_EQ(distance(p1, p2), 0.559016994374947);
  p1 = p2 + 0.5;
  EXPECT_DOUBLE_EQ(p1.x(), 1.0);
  EXPECT_DOUBLE_EQ(p1.y(), 1.0);
  EXPECT_DOUBLE_EQ(p1.z(), 0.5);
  p2 = p1 - 0.1;
  EXPECT_DOUBLE_EQ(p2.x(), 0.9);
  EXPECT_DOUBLE_EQ(p2.y(), 0.9);
  EXPECT_DOUBLE_EQ(p2.z(), 0.4);

 p1 = Point(0.0, 0.0, 0.0);
 p2 = Point(1.0, 0.0, 0.0);
 p3 = Point(1.0, 1.0, 0.0);
 p4 = Point(1.0, 1.0, 1.0);
 auto p5 = 1.0*p1;
 Point p6 {1.0, 0.5, 1.0};
 EXPECT_EQ(p1, p5);
 EXPECT_LT(p1, p2);
 EXPECT_LE(p1, p2);
 EXPECT_LT(p1, p3);
 EXPECT_LE(p1, p3);
 EXPECT_LE(p2, p3);
 EXPECT_LE(p3, p4);
 EXPECT_LT(p6, p4);
}

//---------------------------------------------------------------------------//
//              end of test_Point.cc
//---------------------------------------------------------------------------//
