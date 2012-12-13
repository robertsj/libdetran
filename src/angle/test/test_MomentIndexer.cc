//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_MomentIndexer.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of MomentIndexer class
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_MomentIndexer)

// Detran headers
#include "TestDriver.hh"
#include "MomentIndexer.hh"

// Setup
#include "quadrature_fixture.hh"

using namespace detran_angle;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_MomentIndexer(int argc, char *argv[])
{
  // 1d test
  {
    MomentIndexer indexer(1, 10);
    TEST(indexer.number_moments() == 11);
    TEST(indexer.legendre_order() == 10);
    for (int i = 0; i < 11; ++i)
    {
      TEST(indexer.l(i) == i);
      TEST(indexer.m(i) == 0);
    }
  }
  // 2d test
  {
    int l[] = {0,1,1,2,2,2,3,3,3,3,4,4,4,4,4};
    int m[] = {0,-1,1,-2,0,2,-3,-1,1,3,-4,-2,0,2,4};
    MomentIndexer indexer(2, 4);
    TEST(indexer.number_moments() == 15);
    TEST(indexer.legendre_order() == 4);
    for (int i = 0; i < 15; ++i)
    {
      TEST(indexer.l(i) == l[i]);
      TEST(indexer.m(i) == m[i]);
    }
  }
  // 3d test
  {
    int l[] = {0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4};
    int m[] = {0,-1, 0, 1,-2,-1, 0, 1, 2,-3,-2,-1, 0, 1, 2, 3,-4,-3,-2,-1, 0, 1, 2, 3, 4};
    MomentIndexer indexer(3, 4);
    TEST(indexer.number_moments() == 25);
    TEST(indexer.legendre_order() == 4);
    for (int i = 0; i < 25; ++i)
    {
      TEST(indexer.l(i) == l[i]);
      TEST(indexer.m(i) == m[i]);
    }
  }
  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
