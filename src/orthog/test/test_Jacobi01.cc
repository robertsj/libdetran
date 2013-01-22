//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_Jacobi01.cc
 *  @author Jeremy Roberts
 *  @date   Apr 1, 2012
 *  @brief  Test of Jacobi01 class
 */
//---------------------------------------------------------------------------//

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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

int test_Jacobi01(int argc, char *argv[])
{

  int N = 8;
  int M = 7;

  vec_dbl x(N);
  vec_dbl w(N);
  double xa[] = {-0.9602898564975362316835609,-0.7966664774136267395915539,-0.5255324099163289858177390,-0.1834346424956498049394761,
                  0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609};
  double wa[] = { 0.1012285362903762591525314, 0.2223810344533744705443560, 0.3137066458778872873379622, 0.3626837833783619829651504,
                  0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314};
  for (int i = 0; i < N; ++i)
  {
    x[i] = xa[i];
    w[i] = wa[i];
  }

  Jacobi01 P(M, x, w);
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

  f.display();
  ft.display();
  f2.display();

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_MomentIndexer.cc
//---------------------------------------------------------------------------//
