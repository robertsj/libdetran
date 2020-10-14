//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Preconditioners.cc
 *  @brief Test of Matrix class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST              \
        FUNC(test_PCJacobi)    \
        FUNC(test_PCILU0)      \
        FUNC(test_PCILUT_full) \
        FUNC(test_PCILUT_P0)   \


#include "TestDriver.hh"
#include "preconditioner/PCJacobi.hh"
#include "preconditioner/PCILU0.hh"
#include "preconditioner/PCILUT.hh"

#include "matrix_fixture.hh"
#include "matrix/Matrix.hh"
#include "utils/Initialization.hh"
#include <iostream>

using namespace callow;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
int test_PCJacobi(int argc, char *argv[])
{
  Matrix::SP_matrix A = test_matrix_1(5);

  PCJacobi P(A);
  P.display("pc_jacobi.out");

  return 0;
}

//----------------------------------------------------------------------------//
int test_PCILU0(int argc, char *argv[])
{
  Matrix::SP_matrix A = test_matrix_1(5);

  PCILU0 P(A);
  P.display("pc_ilu0.out");

  return 0;
}


//----------------------------------------------------------------------------//
int test_PCILUT_full(int argc, char *argv[])
{

  {
	// small 2-d diffusion matrix.  By default, ILUT should
	// provide the full LU decomposition.
    Matrix::SP_matrix A = test_matrix_2(2);

    // reference
    const int expect_n = 38;
    const double expect_v[expect_n] = {  2.687837678598e-02,
    		 -5.891883929887e-04, -5.891883929887e-04, -2.192053477337e-02,
    		  2.791961879430e-02, -1.291532465657e-05, -5.891883929887e-04,
    		 -2.192053477337e-02, -4.625895773049e-04,  2.791961281981e-02,
    		 -5.894609453983e-04, -2.110302426869e-02, -2.111279082567e-02,
    		  2.896181262935e-02, -5.989945050699e-01, -1.264059557730e-02,
    		 -1.264644569127e-02, -5.145492175657e-04,  7.666353065191e-02,
    		 -1.317653259545e-04, -1.317653259545e-04, -5.766554378345e-01,
    		 -2.667548523211e-04, -1.173669399231e-02, -1.718748469240e-03,
    		  7.692006831051e-02, -2.264714522831e-07, -1.317653259545e-04,
    		 -5.766555612325e-01, -1.173669399231e-02, -1.718748469240e-03,
    		 -2.944244034845e-06,  7.692006830984e-02, -1.317657139038e-04,
    		 -5.559044320205e-01, -1.713016236837e-03, -1.713021280389e-03,
    		  7.717660747840e-02 };

    PCILUT P(A);
    Matrix::SP_matrix M = P.matrix();

    int got_n = M->number_nonzeros();
    double *got_v = M->values();

    TEST(got_n == expect_n);
    for (int i = 0; i < got_n; ++i)
    {
    	TEST(soft_equiv(got_v[i], expect_v[i]));
    }

    Vector given(M->number_columns(), 1.0);
    Vector got(M->number_columns(), 0.0);
    P.apply(given, got);
    double expect[8] = {3.884339023513e+01, 3.737962161918e+01,
    		            3.737962161918e+01, 3.601815957905e+01,
    		            2.127328717126e+01, 2.089593865027e+01,
						2.089593865027e+01, 2.054236209448e+01};
    for (int i = 0; i < 8; ++i)
    {
    	TEST(soft_equiv(got[i], expect[i]));
    }

  }

  return 0;
}

//----------------------------------------------------------------------------//
int test_PCILUT_P0(int argc, char *argv[])
{

  {
	// small 2-d diffusion matrix.
    Matrix::SP_matrix A = test_matrix_2(2);

    // reference
    const int expect_n = 8;
    const double expect_v[expect_n] = {2.687837678598e-02,	2.793253411896e-02,
			2.793253411896e-02, 2.898669145194e-02,
			7.666353065191e-02, 7.692029478196e-02,
			7.692029478196e-02, 7.717705891201e-02};

    // ILUT(0, 0.0) should yield the diagonal matrix.
    PCILUT P(A, 0, 0.0);
    Matrix::SP_matrix M = P.matrix();
    M->display();
    int got_n = M->number_nonzeros();
    double *got_v = M->values();

    TEST(got_n == expect_n);
    for (int i = 0; i < got_n; ++i)
    {
    	TEST(soft_equiv(got_v[i], expect_v[i]));
    }
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Preconditioners.cc
//----------------------------------------------------------------------------//
