//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Segment.cc
 *  @brief Test of Segment class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//
#include <gtest/gtest.h>
#include "geometry/Segment.hh"

using namespace detran_geometry;
using namespace detran_utilities;

TEST(Segment, Basic)
{
  Segment s(0, 1.0);
  EXPECT_EQ(s.region(), 0);
  EXPECT_NEAR(s.length(), 1.0, 1.0e-12);
}

//----------------------------------------------------------------------------//
//              end of test_Segment.cc
//----------------------------------------------------------------------------//
