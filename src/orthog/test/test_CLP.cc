//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_CLP.cc
 *  @brief Test of CLP class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_CLP)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "utilities/MathUtilities.hh"
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

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_CLP(int argc, char *argv[])
{
  int N = 10;
  int M = 3;

  vec_dbl x = linspace(-1.0, 1.0, N);
  vec_dbl dx(N, x[1]-x[0]);

  CLP P(M, x, dx);
  callow::Vector X(x);
  X.display("X");
  P.basis()->display();

  callow::Vector f(N, 0.0);
  for (int i = 0; i < N; ++i)
    f[i] = x[i]*x[i];

  callow::Vector ft(M + 1, 0.0);
  callow::Vector fa(N,     0.0);

  P.transform(f, ft);

  ft.display("FT");

  P.inverse(ft, fa);

  fa.display("FA");

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_CLP.cc
//----------------------------------------------------------------------------//
