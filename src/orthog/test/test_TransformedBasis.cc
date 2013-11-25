//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_TransformedBasis.cc
 *  @brief Test of TransformedBasis class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_TransformedBasis)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "orthog/TransformedBasis.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>

using namespace detran_orthog;
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

int test_TransformedBasis(int argc, char *argv[])
{
  callow_initialize(argc, argv);

  {
    TransformedBasis::Parameters p;
    p.size = 10;
    p.order = 3;
    p.orthonormal = true;
    p.x.resize(p.size, 0.0);
    double x[] = {7, 10, 6, 1, 0.5, 0.2, 0.1, 0.05, 0.001, 0.0001};
    for (int i = 0; i < 10; ++i)
      p.x[i] = x[i];

    OrthogonalBasis::Factory_T::ShowRegistered();
    OrthogonalBasis::SP_basis P = OrthogonalBasis::Create("trans", p);
    P->basis()->print_matlab("P.out");

    double ref[] = {5.128480369198469e-01, 1.315046629466771e-01,
        3.230793420994966e-01, 3.896332471932851e-01};

    for (int i = 0; i < 4; ++i)
    {
      TEST(soft_equiv((*P->basis())(0, i), ref[i]));
    }

    double ref_coef[] = {6.432239496678253e-01, -1.967195237976410e-01,
        5.457188438895422e-03, 2.181167781209342e+00};

    double ref_appx[] = {1.155625195331455e+00, 4.848937080558042e-01,
        -5.211589957765328e-01, -9.819183888198760e-01, -5.030854757398791e-01,
        3.768960079559328e-01, 9.220275148302088e-01, 6.852283179257563e-01,
        -1.239559670944787e-01, -8.206175147406938e-01};

    callow::Vector f(p.size, 0.0);
    for (int i = 0; i < p.size; ++i)
      f[i] = std::cos(i);
    callow::Vector ft(p.order + 1, 0.0);
    callow::Vector fa(p.size, 0.0);

    P->transform(f, ft);

    f.display("F");
    ft.display("FT");
    for (int i = 0; i < 4; ++i)
    {
      TEST(soft_equiv(ft[i], ref_coef[i]));
    }

    P->inverse(ft, fa);

    fa.display("FA");
    for (int i = 0; i < 10; ++i)
    {
      TEST(soft_equiv(fa[i], ref_appx[i]));
    }
  }

  callow_finalize();

  return 0;
}


//----------------------------------------------------------------------------//
//              end of test_TransformedBasis.cc
//----------------------------------------------------------------------------//
