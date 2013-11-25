//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_PolarQuadrature.cc
 *  @brief Test of GaussLegendre class
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_PolarQuadrature)

#include "TestDriver.hh"
#include "QuadratureFactory.hh"

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

// Make sure the a few quadratures give expected numbers
int test_PolarQuadrature(int argc, char *argv[])
{
  int N = 4;
  const char* types[] = {"gl", "dgl", "gc", "gc"};
  double mu[4][3] = {{0.1252334085115, 0.3678314989982, 0.5873179542866},
                     {0.0337652428984, 0.1693953067669, 0.3806904069584},
                     {0.1305261922201, 0.3826834323651, 0.6087614290087},
                     {0.1305261922201, 0.3826834323651, 0.6087614290087}};
  double wt[4][3] = {{0.2491470458134, 0.2334925365384, 0.2031674267231},
                     {0.0856622461896, 0.1803807865241, 0.2339569672863},
                     {0.2595596577443, 0.2418710960116, 0.2076994187967},
                     {0.2588190451025, 0.2411809548975, 0.2071067811865}};
  int norm[4] = {0, 0, 0, 1};

  InputDB::SP_input db = InputDB::Create();
  db->put<int>("quad_number_polar_octant",   6);
  for (int i = 0; i < N; ++i)
  {
    db->put<std::string>("quad_type", types[i]);
    db->put<int>("quad_chebyshev_normalize", norm[i]);
    QuadratureFactory::SP_quadrature q = QuadratureFactory::build(db, 1);
    for (int j = 0; j < 3; ++j)
    {
      TEST(soft_equiv(mu[i][j], q->mu(0, j),  1e-11));
      TEST(soft_equiv(wt[i][j], q->weight(j), 1e-11));
    }
  }

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_PolarQuadrature.cc
//---------------------------------------------------------------------------//
