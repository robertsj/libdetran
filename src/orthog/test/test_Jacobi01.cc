//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Jacobi01.cc
 *  @brief Test of Jacobi01 class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_Jacobi01)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/Jacobi01.hh"
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

int test_Jacobi01(int argc, char *argv[])
{

  int N = 4;
  int M = 3;

  vec_dbl x(N);
  vec_dbl w(N);
  double xa[] = {0.9305681557970262, 0.6699905217924281, 0.33000947820757187, 0.06943184420297371};
  double wa[] = {0.17392742256872687, 0.3260725774312732, 0.3260725774312732, 0.17392742256872687};
  for (int i = 0; i < N; ++i)
  {
    x[i] = xa[i];
    w[i] = wa[i];
  }

  Jacobi01 P(M, x, w, 0.0, 1.0);
  P.basis()->print_matlab("P.out");
  P.weights()->print_matlab("W.out");
  P.coefficients()->print_matlab("A.out");

  callow::Vector f(N, 0.0);
  for (int i = 0; i < N; ++i) f[i] = std::cos(x[i]);
  f.print_matlab("F.out");
  callow::Vector ft(M + 1, 0.0);
  callow::Vector f2(N, 0.0);

  P.transform(f, ft);
  P.inverse(ft, f2);

  std::cout << " **** " << std::endl;
  f.display();
  ft.display();
  f2.display();
  std::cout << " L2 res norm = " << f2.norm_residual(f) << std::endl;

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_Jacobi01.cc
//----------------------------------------------------------------------------//
