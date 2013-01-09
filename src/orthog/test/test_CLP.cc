//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_CLP.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of CLP class
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_CLP)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/CLP.hh"
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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_CLP(int argc, char *argv[])
{

  int N = 20;
  int M = 4;
  double width = 2.0 / N;

  vec_dbl x(N, 0);
  x[0] = -1.0 + width / 2.0;
  for (int i = 1; i < N; ++i) x[i] = x[i-1] + width;

  vec_dbl dx(N, width);

  CLP P(M, x, dx);

  callow::Vector f(N, 0.0);
  for (int i = 0; i < N; ++i) f[i] = std::cos(x[i]);

  callow::Vector ft(M + 1, 0.0);
  callow::Vector f2(N, 0.0);

  P.transform(f, ft);
  P.inverse(ft, f2);

  f.display();
  ft.display();
  f2.display();


  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
