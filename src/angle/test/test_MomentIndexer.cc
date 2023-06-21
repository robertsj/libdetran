//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MomentIndexer.cc
 *  @brief Test of MomentIndexer class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "MomentIndexer.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace std;

TEST(MomentIndexer, 1D)
{
  std::cout << "1D case" << std::endl;
  MomentIndexer indexer(1, 10);
  EXPECT_EQ(indexer.number_moments(), 11);
  EXPECT_EQ(indexer.legendre_order(), 10);
  for (int i = 0; i < 11; ++i)
  {
    EXPECT_EQ(indexer.l(i), i);
    EXPECT_EQ(indexer.m(i), 0);
    EXPECT_EQ(indexer.index(i, 0), i);
  }

}

TEST(MomentIndexer, 2D)
{
  int l[] =
  { 0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4 };
  int m[] =
  { 0, -1, 1, -2, 0, 2, -3, -1, 1, 3, -4, -2, 0, 2, 4 };
  MomentIndexer indexer(2, 4);
  EXPECT_EQ(indexer.number_moments(), 15);
  EXPECT_EQ(indexer.legendre_order(), 4);
  for (int i = 0; i < 15; ++i)
  {
    EXPECT_EQ(indexer.l(i), l[i]);
    EXPECT_EQ(indexer.m(i), m[i]);
    EXPECT_EQ(indexer.index(l[i], m[i]), i);
  }

}

TEST(MomentIndexer, 3D)
{
  int l[] =
  { 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4 };
  int m[] =
  { 0, -1, 0, 1, -2, -1, 0, 1, 2, -3, -2, -1, 0, 1, 2, 3, -4, -3, -2, -1, 0, 1,
      2, 3, 4 };
  MomentIndexer indexer(3, 4);
  EXPECT_EQ(indexer.number_moments(), 25);
  EXPECT_EQ(indexer.legendre_order(), 4);
  for (int i = 0; i < 25; ++i)
  {
    EXPECT_EQ(indexer.l(i), l[i]);
    EXPECT_EQ(indexer.m(i), m[i]);
    EXPECT_EQ(indexer.index(l[i], m[i]), i);
  }

}

//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
