//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_ChebyshevU.cc
 *  @brief Test of ChebyshevU class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_ChebyshevU)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
#include "utilities/MathUtilities.hh"
#include "orthog/ChebyshevU.hh"
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

int test_ChebyshevU(int argc, char *argv[])
{

  {
    ChebyshevU::Parameters p;
    p.size = 10;
    p.order = 6;
    p.x = linspace(-1, 1.0, (int)p.size);
    p.qw = vec_dbl(p.size, p.x[1]-p.x[0]);

    ChebyshevU P(p);

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
  }

  {
    ChebyshevU::Parameters p;
    p.size = 4;
    p.order = 3;
    p.x.resize(p.size, 0.0);
    p.qw.resize(p.size, 0.0);
    p.lower_bound = 0.0;
    p.upper_bound = 1.0;

    double xa[] = {0.9619397662556434, 0.6913417161825449, 0.30865828381745514, 0.03806023374435663};
    double wa[] = {0.15027943247108658, 0.36280664401742885, 0.36280664401742885, 0.15027943247108658};
//    double xa[] = {0.9305681557970262, 0.6699905217924281, 0.33000947820757187, 0.06943184420297371};
//    double wa[] = {0.17392742256872687, 0.3260725774312732, 0.3260725774312732, 0.17392742256872687};
    for (int i = 0; i < p.size; ++i)
    {
      p.x[i] = xa[i];
      p.qw[i] = wa[i];
    }

    ChebyshevU P(p);

    callow::Vector f(p.size, 0.0);
    for (int i = 0; i < p.size; ++i)
      f[i] = std::cos(p.x[i]);
    callow::Vector ft(p.order + 1, 0.0);
    callow::Vector fa(p.size, 0.0);

    P.transform(f, ft);
    P.inverse(ft, fa);

    std::cout << " **** " << std::endl;
    f.display("F");
    ft.display("FT");
    fa.display("FA");
    std::cout << " L2 res norm = " << fa.norm_residual(f) << std::endl;

    TEST(soft_equiv(fa.norm_residual(f), 0.0));
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_ChebyshevU.cc
//----------------------------------------------------------------------------//
