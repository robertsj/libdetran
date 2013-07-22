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
  CLP::Parameters p;
  p.size = 10;
  p.order = 3;
  p.orthonormal = false;
  p.x = linspace_center(-1.0, 1.0, p.size);
  p.qw = vec_dbl(p.size, p.x[1]-p.x[0]);

  callow::Vector X(p.x);
  X.display("X");

  OrthogonalBasis::Factory_T::ShowRegistered();
  OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("clp", p);
  P->basis()->display();

  callow::Vector f(p.size, 0.0);
  for (int i = 0; i < p.size; ++i)
    f[i] = p.x[i]*p.x[i];

  callow::Vector ft(p.order+ 1, 0.0);
  callow::Vector fa(p.size,     0.0);

  double ref_coef[] = { 6.6e-01, 0.0, 2.5014e-01, 0.0};

  P->transform(f, ft);

  f.display("F");
  ft.display("FT");

  for (int i = 0; i <= p.order; ++i)
  {
    TEST(soft_equiv(ft[i], ref_coef[i]));
  }

  double ref_appx[] = {7.771252499999999e-01, 4.769572499999999e-01,
      2.518312500000000e-01, 1.017472500000000e-01, 2.670525000000006e-02,
      2.670525000000008e-02, 1.017472500000000e-01, 2.518312500000002e-01,
      4.769572500000001e-01, 7.771252500000000e-01};

  P->inverse(ft, fa);

  fa.display("FA");

  for (int i = 0; i < p.size; ++i)
  {
    TEST(soft_equiv(fa[i], ref_appx[i]));
  }

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_CLP.cc
//----------------------------------------------------------------------------//
