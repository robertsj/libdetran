//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_DCT.cc
 *  @brief Test of DCT class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_DCT)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "utilities/MathUtilities.hh"
#include "orthog/DCT.hh"
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

int test_DCT(int argc, char *argv[])
{
  DCT::Parameters p;
  p.size  = 10;
  p.order = 3;

  // orthonormal
  {
    p.orthonormal = true;

    DCT P(p);

    callow::Vector f(p.size, 0.0);
    for (int i = 0; i < p.size; ++i)
      f[i] = std::cos(i);

    callow::Vector ft(p.order + 1, 0.0);
    callow::Vector fa(p.size, 0.0);

    P.basis()->display();

    P.transform(f, ft);

    double ref_coef[] = {1.333291468776643e-01, 2.859141524263383e-01,
        3.067802114121704e-01, 2.231212702582738e+00};

    P.inverse(ft, fa);

    double ref_appx[] = {1.188006106865787e+00, 3.928273994252888e-01,
        -5.729950364896974e-01, -9.659739438142930e-01, -5.213213235757106e-01,
        3.446832452756474e-01, 8.890147986165391e-01, 6.573197930138066e-01,
        -1.472187411793163e-01, -8.427185155175049e-01};

    for (int i = 0; i < 4; ++i)
    {
      TEST(soft_equiv(ft[i], ref_coef[i]));
    }

    for (int i = 0; i < 10; ++i)
    {
      TEST(soft_equiv(fa[i], ref_appx[i]));
    }
  }


  // not orthonormal
  {
    p.orthonormal = false;

    DCT P(p);
    P.coefficients()->display("A");

    callow::Vector f(p.size, 0.0);
    for (int i = 0; i < p.size; ++i)
      f[i] = std::cos(i);

    callow::Vector ft(p.order + 1, 0.0);
    callow::Vector fa(p.size, 0.0);

    P.basis()->display();

    P.transform(f, ft);

    double ref_coef[] = {4.216237826205462e-01, 6.393234805545289e-01,
        6.859814068693698e-01, 4.989143275236023e+00};

    P.inverse(ft, fa);

    double ref_appx[] = {1.188006106865787e+00, 3.928273994252888e-01,
        -5.729950364896974e-01, -9.659739438142930e-01, -5.213213235757106e-01,
        3.446832452756474e-01, 8.890147986165391e-01, 6.573197930138066e-01,
        -1.472187411793163e-01, -8.427185155175049e-01};

    for (int i = 0; i < 4; ++i)
    {
      TEST(soft_equiv(ft[i], ref_coef[i]));
    }

    for (int i = 0; i < 10; ++i)
    {
      TEST(soft_equiv(fa[i], ref_appx[i]));
    }
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_DCT.cc
//----------------------------------------------------------------------------//
