//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MomentIndexer.cc
 *  @brief Test of MomentIndexer class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_MomentIndexer_1D) \
        FUNC(test_MomentIndexer_2D) \
        FUNC(test_MomentIndexer_3D)

#include "TestDriver.hh"
#include "MomentIndexer.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_MomentIndexer_1D(int argc, char *argv[])
{
  std::cout << "1D case" << std::endl;
  MomentIndexer indexer(1, 10);
  TEST(indexer.number_moments() == 11);
  TEST(indexer.legendre_order() == 10);
  for (int i = 0; i < 11; ++i)
  {
    TEST(indexer.l(i) == i);
    TEST(indexer.m(i) == 0);
    TEST(indexer.index(i, 0) == i);
  }

  return 0;
}

int test_MomentIndexer_2D(int argc, char *argv[])
{
  int l[] =
  { 0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4 };
  int m[] =
  { 0, -1, 1, -2, 0, 2, -3, -1, 1, 3, -4, -2, 0, 2, 4 };
  MomentIndexer indexer(2, 4);
  TEST(indexer.number_moments() == 15);
  TEST(indexer.legendre_order() == 4);
  for (int i = 0; i < 15; ++i)
  {
    TEST(indexer.l(i) == l[i]);
    TEST(indexer.m(i) == m[i]);
    TEST(indexer.index(l[i], m[i]) == i);
  }

  return 0;
}

int test_MomentIndexer_3D(int argc, char *argv[])
{
  int l[] =
  { 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4 };
  int m[] =
  { 0, -1, 0, 1, -2, -1, 0, 1, 2, -3, -2, -1, 0, 1, 2, 3, -4, -3, -2, -1, 0, 1,
      2, 3, 4 };
  MomentIndexer indexer(3, 4);
  TEST(indexer.number_moments() == 25);
  TEST(indexer.legendre_order() == 4);
  for (int i = 0; i < 25; ++i)
  {
    TEST(indexer.l(i) == l[i]);
    TEST(indexer.m(i) == m[i]);
    TEST(indexer.index(l[i], m[i]) == i);
  }

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
