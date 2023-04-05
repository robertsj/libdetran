//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_TinyVector.cc
 *  @brief Test of InputDB
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "utilities/TinyVector.hh"

using namespace detran_utilities;

TEST(TinyVector, Basic)
{
  TinyVector<int, 2> V2(1, 2);
  EXPECT_EQ(V2.size(), 2);
  EXPECT_EQ(V2[0], 1);
  EXPECT_EQ(V2[1], 2);

  TinyVector<int, 3> V3(1, 2, 3);
  EXPECT_EQ(V3.size(), 3);
  EXPECT_EQ(V3[0], 1);
  EXPECT_EQ(V3[1], 2);
  EXPECT_EQ(V3[2], 3);
}

//----------------------------------------------------------------------------//
//              end of test_TinyVector.cc
//----------------------------------------------------------------------------//
