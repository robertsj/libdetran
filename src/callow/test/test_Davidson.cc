//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Davidson.cc
 *  @brief Test of test_Davidson class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "callow/utils/Initialization.hh"
#include "callow/solver/Davidson.hh"
#include "callow/solver/PowerIteration.hh"

#include "matrix_fixture.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

class DiagonalPC: public PCShell
{
public:
  DiagonalPC(Matrix::SP_matrix A, Matrix::SP_matrix B, void* context)
    : PCShell("diag", context), d_A(A), d_B(B)
  {/* ... */}
  void apply(Vector &x, Vector &y)
  {
    std::cout << " applied" << std::endl;
    double e = ((Davidson*)PCShell::d_context)->eigenvalue();
    for (size_t i = 0; i < x.size(); ++i)
    {
      double val = (*d_A)[d_A->diagonal(i)] - 1 * (*d_B)[d_B->diagonal(i)];
      y[i] = x[i] / val;
    }
  }
private:
  Matrix::SP_matrix d_A;
  Matrix::SP_matrix d_B;
};

/*
 *  Solves Ax = e*B*x where A is the standard 2nd difference
 *  operator and B is diagonal matrix of 1 through n.
 */
TEST(Davidson, Standard)
{
  Matrix::SP_matrix M = test_matrix_2(10);
  Matrix::SP_matrix F = test_matrix_3(10);

  M->print_matlab("M.out");
  F->print_matlab("F.out");

  int n = M->number_columns();
  Vector X (n, 1.0);
  Vector X0(n, 1.0);

  // test on Ax=Lx
  {
    double ref_L = 0.141079562648794;
    double ref_V[] =
    { 1.000000000000000, -2.916559254113351, 4.589758628640770};
    Davidson solver(1e-10, 100, 20);
    solver.set_operators(M);
    solver.solve(X, X0);
    X.scale(1.0/X[0]);
    printf("%16.8f  %16.8f \n", solver.eigenvalue(), ref_L);
    EXPECT_NEAR(solver.eigenvalue(), ref_L, 1.0e-8);
    for (int i = 0; i < 3; ++i)
    {
      printf("%16.8f  %16.8f \n", X[i], ref_V[i]);
      EXPECT_NEAR(X[i], ref_V[i], 1.0e-8);
    }
  }  
}

TEST(Davidson, General)
{
  Matrix::SP_matrix M = test_matrix_2(10);
  Matrix::SP_matrix F = test_matrix_3(10);
  int n = M->number_columns();
  Vector X (n, 1.0);
  Vector X0(n, 1.0);

  // Test on Ax=LBx
  {
    double ref_L = 1.243023126562274;
    double ref_V[] =
    { 1.000000000000000, 0.976684255089969, 0.930596388729471};
    Davidson solver(1e-10, 50, 20);
    solver.set_operators(F, M);
    solver.solve(X, X0);
    X.scale(1.0/X[0]);
    printf("%16.8f  %16.8f \n", solver.eigenvalue(), ref_L);
    EXPECT_NEAR(solver.eigenvalue(), ref_L, 1.0e-8);
    for (int i = 0; i < 3; ++i)
    {
      printf("%16.8f  %16.8f \n", X[i], ref_V[i]);
      EXPECT_NEAR(X[i], ref_V[i], 1.0e-8);
    }
  }
}

//---------------------------------------------------------------------------//
//              end of test_Matrix.cc
//---------------------------------------------------------------------------//
