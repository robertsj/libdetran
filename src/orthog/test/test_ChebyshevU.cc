//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_ChebyshevU.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of CLP class
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                   \
        FUNC(test_ChebyshevU)

#include "utilities/TestDriver.hh"
#include "utilities/Definitions.hh"
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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_ChebyshevU(int argc, char *argv[])
{

  if(0){
    int N = 10;
    int M = 6;
    double width = 2.0 / N;

    vec_dbl x(N, 0);
    x[0] = -1.0 + width / 2.0;
    for (int i = 1; i < N; ++i) x[i] = x[i-1] + width;

    vec_dbl dx(N, width);

    ChebyshevU P(M, x, dx);

    callow::Vector f(N, 0.0);
    for (int i = 0; i < N; ++i) f[i] = std::cos(x[i]);

    callow::Vector ft(M + 1, 0.0);
    callow::Vector f2(N, 0.0);

    P.transform(f, ft);
    P.inverse(ft, f2);

    f.display();
    ft.display();
    f2.display();
  }

  {
    int N = 4;
    int M = 3;

    vec_dbl x(N);
    vec_dbl w(N);
    double xa[] = {0.9619397662556434, 0.6913417161825449, 0.30865828381745514, 0.03806023374435663};
    double wa[] = {0.15027943247108658, 0.36280664401742885, 0.36280664401742885, 0.15027943247108658};
//    double xa[] = {0.9305681557970262, 0.6699905217924281, 0.33000947820757187, 0.06943184420297371};
//    double wa[] = {0.17392742256872687, 0.3260725774312732, 0.3260725774312732, 0.17392742256872687};
    for (int i = 0; i < N; ++i)
    {
      x[i] = xa[i];
      w[i] = wa[i];
    }

    ChebyshevU P(M, x, w, 0.0, 1.0);

    callow::Vector f(N, 0.0);
    for (int i = 0; i < N; ++i) f[i] = std::cos(x[i]);
    callow::Vector ft(M + 1, 0.0);
    callow::Vector f2(N, 0.0);

    P.transform(f, ft);
    P.inverse(ft, f2);

    std::cout << " **** " << std::endl;
    f.display();
    ft.display();
    f2.display();
    std::cout << " L2 res norm = " << f2.norm_residual(f) << std::endl;

  }

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
