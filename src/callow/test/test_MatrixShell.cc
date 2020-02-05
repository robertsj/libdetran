//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_MatrixShell.cc
 *  @brief Test of MatrixShell class
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST              \
        FUNC(test_MatrixShell)

#include "TestDriver.hh"
#include "matrix/MatrixShell.hh"
#include "utils/Initialization.hh"
#include "test/matrixshell_fixture.hh"
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

// Test of basic public interface
int test_MatrixShell(int argc, char *argv[])
{
  typedef TestMatrixShell Mat_T;
  typedef Vector          Vec_T;

  // n * n
  {
    // Create test matrix
    int n = 5;
    Mat_T A(n);
    A.display();
    TEST(A.number_rows()    == n);
    TEST(A.number_columns() == n);

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
      TEST(soft_equiv(Y[i], ref[i]));
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    {-0.037,  -0.1079,  -0.0416, -1.0511,   3.0286};

    for (int i = 0; i < X.size(); i++)
    {
      TEST(soft_equiv(X[i], ref2[i]));
    }

  } // end n * n

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_MatrixShell.cc
//----------------------------------------------------------------------------//
