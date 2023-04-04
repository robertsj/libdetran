//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Matrix.cc
 *  @brief Test of Matrix class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "matrix_fixture.hh"
#include "matrix/Matrix.hh"
#include "utils/Initialization.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

// Test of basic public interface
TEST(Matrix, Basic)
{
  // n * n
  {
    // Create test matrix
    int n = 5;
    Matrix A(n, n);
    EXPECT_EQ(A.number_rows()   , n);
    EXPECT_EQ(A.number_columns(), n);
    A.preallocate(3);

    for (int row = n - 1; row >= 0; row--)
    {
      if (row == 0)
      {
        int c[]    = { 0, 1  };
        double v[] = {-2, 1.1};
        EXPECT_TRUE(A.insert(row, c, v, 2));
      }
      else if (row == A.number_rows() - 1)
      {
        int c[] = {A.number_rows() - 2, A.number_rows() - 1};
        double v[] = {1, -2};
        EXPECT_TRUE(A.insert(row, c, v, 2));
      }
      else
      {
        int c[] = {row - 1, row, row + 1};
        double v[] = {1, -2, 1.1};
        EXPECT_TRUE(A.insert(row, c, v, 3));
      }
    }
    A.assemble();
    A.display();
    // Ensure no repeated columns
    for (int i = 0; i < A.number_rows(); ++i)
    {
      int c0 = -1, c1 = 0;
      for (int p = A.start(i); p < A.end(i); ++p)
      {
        c1 = A.column(p);
        EXPECT_TRUE(c1 > c0);
        c0 = c1;
      }
    }

    // Create two vectors
    Vector X(n, 0.0);
    Vector Y(n, 0.0);

    for (int i = 0; i < X.size(); i++)
    {
      double value = i * i;
      X[i] = value;
    }

    // Multiply
    A.multiply(X, Y);
    double ref[] =
    { 1.1, 2.4, 2.9, 3.6, -23.0 };
    for (int i = 0; i < Y.size(); i++)
    {
      EXPECT_NEAR(Y[i], ref[i], 1.0e-12);
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    { 0.2, -0.69, 0.44, -27.01, 49.96 };
    for (int i = 0; i < X.size(); i++)
    {
      EXPECT_NEAR(X[i], ref2[i], 1.0e-12);
    }

    // Indexing
    for (int i = 0; i < A.number_rows(); i++)
    {
      std::cout << i << " " << A.diagonal(i) << std::endl;
//      for (int p = A.start(i); p < A.diagonal(i); p++)
//      {
//        std::cout << " " << A.column(p) << "(" << a[p] << ")";
//      }
//      std::cout << endl;
    }

    A.display();
    X.display();
    Y.display();

  } // end n * n

  // m * n
  {
    int m = 3;
    int n = 2;

    // Create test matrix
    Matrix A(m, n);
    EXPECT_EQ(A.number_rows()   , m);
    EXPECT_EQ(A.number_columns(), n);
    A.preallocate(2);
    /*!
     *  1 2
     *  3 4
     *  5 6
     */
    A.insert(0, 0, 1.0);
    A.insert(0, 1, 2.0);
    A.insert(1, 0, 3.0);
    A.insert(1, 1, 4.0);
    A.insert(2, 0, 5.0);
    A.insert(2, 1, 6.0);
    A.assemble();

    // Create two vectors
    Vector X(n, 1.0);
    Vector Y(m, 0.0);

    // Multiply
    A.multiply(X, Y);
    double ref[] =
    { 3.0, 7.0, 11.0 };
    for (int i = 0; i < Y.size(); i++)
    {
      EXPECT_NEAR(Y[i], ref[i], 1.0e-12);
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    { 79.0, 100.0 };
    for (int i = 0; i < X.size(); i++)
    {
      EXPECT_NEAR(X[i], ref2[i], 1.0e-12);
    }

    A.display();
    X.display();
    Y.display();

  } // end m * n

  // Pathologic case?
  {
    std::cout << " patho case 1" << std::endl;
    // | 1  0 |
    // | 1  0 |
    // | 0  1 |
    // | 0  1 |
    Matrix A(4, 2, 1);
    A.insert(0, 0, 1.0); A.insert(1, 0, 1.0);
    A.insert(2, 1, 1.0); A.insert(3, 1, 1.0);
    A.assemble();
    A.display();
  }

  // Pathologic case?
  {
    std::cout << " patho case 2" << std::endl;
    // | 1  1  0  1 |
    // | 0  0  1  1 |
    Matrix A(2, 4, 2);
    A.insert(0, 0, 1.0); A.insert(0, 1, 1.0);
    A.insert(1, 2, 1.0); A.insert(1, 3, 1.0);
    A.assemble();
    A.display();
  }

  // Clearing
  {
    Matrix A(3, 3, 2);
    A.insert(0, 0, 1.0); A.insert(0, 1, 2.0);
    A.insert(1, 1, 3.0); A.insert(1, 2, 4.0);
    A.insert(2, 1, 5.0); A.insert(2, 2, 6.0);
    A.assemble();
    A.display();
    A.clear();
    A.insert(0, 1, 1.1); A.insert(0, 2, 2.1);
    A.insert(1, 0, 3.1); A.insert(1, 2, 4.1);
    A.insert(2, 0, 5.1); A.insert(2, 2, 6.1);
    A.assemble();
    A.display();
  }
  callow_finalize();
}

TEST(Matrix, Diff)
{
  Matrix::SP_matrix L = test_matrix_2(10);
  Matrix::SP_matrix F = test_matrix_3(10);
  L->print_matlab("L.out");
  F->print_matlab("F.out");
}

//----------------------------------------------------------------------------//
//              end of test_Matrix.cc
//----------------------------------------------------------------------------//
