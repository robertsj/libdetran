//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_DCP.cc
 *  @brief Test of DCP class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_DCP)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/DCP.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_orthog;
using namespace detran_utilities;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_DCP(int argc, char *argv[])
{
  int N = 10;   // size
  int M = 3;   // order
  double width = 2.0 / N;

  DCP P(M, N);

  double ref[] = {0.316227766016838, 0.495433694306862, 0.522232967867093,
      0.453425192941883};

  for (int i = 0; i < 4; ++i)
  {
   // TEST(soft_equiv((*P.basis())(0, i), ref[i]));
  }

  double ref_coef[] = {0.133329146877664, 0.515402376383869, 0.222337198116035,
      1.979100267734324};

  double ref_appx[] = {1.310995817236610, -0.019654554407243,
      -0.583141356287192, -0.593125045724381, -0.263266080039961,
      0.192775083444925, 0.561337987409130, 0.628762174531506,
      0.181387187490909, -0.994447431033812};

  callow::Vector f(N, 0.0);
  for (int i = 0; i < N; ++i)
    f[i] = std::cos(i);
  callow::Vector ft(M + 1, 0.0);
  callow::Vector fa(N, 0.0);

  P.transform(f, ft);

  ft.display("FT");
  for (int i = 0; i < 4; ++i)
  {
   // TEST(soft_equiv(ft[i], ref_coef[i]));
  }


  P.inverse(ft, fa);

  fa.display("FA");
  for (int i = 0; i < 10; ++i)
  {
   // TEST(soft_equiv(fa[i], ref_appx[i]));
  }

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_DCP.cc
//----------------------------------------------------------------------------//
