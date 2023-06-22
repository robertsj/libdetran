//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MatrixShell.cc
 *  @brief Test of MatrixShell class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "matrix/MatrixShell.hh"
#include "utils/Initialization.hh"
#include "test/matrixshell_fixture.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
TEST(MatrixShell, Basic)
{
  typedef TestMatrixShell Mat_T;
  typedef Vector          Vec_T;

  // n * n
  {
    // Create test matrix
    int n = 5;
    Mat_T A(n);
    A.display();
    EXPECT_EQ(A.number_rows()   , n);
    EXPECT_EQ(A.number_columns(), n);

    // Create two vectors
    Vec_T X(n, 0.0);
    Vec_T Y(n, 0.0);

    for (int i = 0; i < X.size(); i++)
    {
      double value = i * i;
      X[i] = value;
    }

    // Multiply
    A.multiply(X, Y);
    double ref[] =
    { -0.21,  -0.34,  -0.09,   0.34,   6.2};
    Y.display();
    for (int i = 0; i < Y.size(); i++)
    {
      EXPECT_NEAR(Y[i], ref[i], 1.0e-12);
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    {-0.037,  -0.1079,  -0.0416, -1.0511,   3.0286};

    for (int i = 0; i < X.size(); i++)
    {
      EXPECT_NEAR(X[i], ref2[i], 1.0e-12);
    }

  } // end n * n
}

// Test matlab writer
TEST(MatrixShell, Write)
{
  // n * n
  {
    // Create test matrix
    int n = 5;
    TestMatrixShell A(n);
    A.print_matlab("shell_matlab.out");
  } // end n * n
}

//----------------------------------------------------------------------------//
//              end of test_MatrixShell.cc
//----------------------------------------------------------------------------//
