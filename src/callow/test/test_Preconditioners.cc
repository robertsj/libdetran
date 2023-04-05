//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Preconditioners.cc
 *  @brief Test of Matrix class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "preconditioner/PCJacobi.hh"
#include "preconditioner/PCILU0.hh"
#include "preconditioner/PCILUT.hh"

#include "matrix_fixture.hh"
#include "matrix/Matrix.hh"
#include "solver/GMRES.hh"
#include "utils/Initialization.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
TEST(Preconditioner, PCJacobi)
{
  Matrix::SP_matrix A = test_matrix_1(5);
  PCJacobi P(A);
  P.display("pc_jacobi.out");
}

//----------------------------------------------------------------------------//
TEST(Preconditioner, PCILU0)
{
  Matrix::SP_matrix A = test_matrix_1(5);
  PCILU0 P(A);
  P.display("pc_ilu0.out");
}

//----------------------------------------------------------------------------//
TEST(Preconditioner, FullPCILUT)
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

    PCILUT P(A, 1000, 0.0);
    Matrix::SP_matrix M = P.matrix();
	M->print_matlab("ilu_oo.out");

    int got_n = M->number_nonzeros();
    double *got_v = M->values();

    EXPECT_EQ(got_n, expect_n);
    for (int i = 0; i < got_n; ++i)
    {
    	EXPECT_NEAR(got_v[i], expect_v[i], 1.0e-12);
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
    	EXPECT_NEAR(got[i], expect[i], 1.0e-12);
    }
}

//----------------------------------------------------------------------------//
TEST(Preconditioner, P0PCILUT)
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

    EXPECT_EQ(got_n, expect_n);
    for (int i = 0; i < got_n; ++i)
    {
    	EXPECT_NEAR(got_v[i], expect_v[i], 1.0e-12);
	}    
}

//----------------------------------------------------------------------------//
TEST(Preconditioner, Performance)
{
	/*
	 *  Test the performance of the difference built-in PC's with
	 *  GMRES to ensure things "seem" right.
	 */
	Matrix::SP_matrix A = test_matrix_2(20);
	Vector b(A->number_rows(), 1.0);
	Vector x(A->number_rows(), 0.0);
	GMRES solver(1e-12, 1e-12, 1000, 20);

	A->print_matlab("A_pc.out");

	solver.set_operators(A);

	// No PC
	solver.solve(b, x);

	// PC ILU0
	{
	x.set(0.0);
	auto P = std::make_shared<PCILU0>(A);
	P->matrix()->print_matlab("ilu0.out");
	solver.set_preconditioner(P, 1);
	solver.solve(b, x);
	}

	// PC ILUT(oo, 0)
	{
	x.set(0.0);
	auto P = std::make_shared<PCILUT>(A, 10000, 0.0);
	P->matrix()->print_matlab("ilutoo.out");
	solver.set_preconditioner(P, 1);
	solver.solve(b, x);
	}

	// PC ILUT(0, 0)
	{
	x.set(0.0);
	auto P = std::make_shared<PCILUT>(A, 0, 0.0);
	P->matrix()->print_matlab("ilut0.out");
	solver.set_preconditioner(P, 1);
	solver.solve(b, x);
	}

	// PC ILUT(10, 1e-3)
	{
	x.set(0.0);
	auto P = std::make_shared<PCILUT>(A, 10, 1e-3);
	P->matrix()->print_matlab("ilut2.out");
	solver.set_preconditioner(P, 1);
	solver.solve(b, x);
	}	
}

//----------------------------------------------------------------------------//
//              end of test_Preconditioners.cc
//----------------------------------------------------------------------------//
