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
  // This tests expansion on something Jacobi should get exactly right

  Jacobi01::Parameters p;
  p.size = 4;
  p.order = 3;
  p.orthonormal = false;
  p.x.resize(p.size, 0.0);
  p.qw.resize(p.size, 0.0);
  p.lower_bound = 0.0;
  p.upper_bound = 1.0;
  double xa[] = {0.9305681557970262, 0.6699905217924281, 0.33000947820757187,
      0.06943184420297371};
  double wa[] = {0.17392742256872687, 0.3260725774312732, 0.3260725774312732,
      0.17392742256872687};
  for (int i = 0; i < p.size; ++i)
  {
    p.x[i] = xa[i];
    p.qw[i] = wa[i];
  }

  Jacobi01 P(p);

  callow::Vector f(p.size, 0.0);
  for (int i = 0; i < p.size; ++i)
    f[i] = std::cos(p.x[i]);

  callow::Vector ft(p.order + 1, 0.0);
  callow::Vector fa(p.size, 0.0);

  P.transform(f, ft);
  P.inverse(ft, fa);

  f.display("F");
  ft.display("FT");
  fa.display("FA");

  double fa_norm = fa.norm_residual(f);

  TEST(soft_equiv(fa_norm, 0.0));

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_Jacobi01.cc
//----------------------------------------------------------------------------//
