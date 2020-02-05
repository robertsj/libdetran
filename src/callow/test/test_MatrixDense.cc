//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   test_MatrixDense.cc
 *  @author Jeremy Roberts
 *  @date   Aug 19, 2012
 *  @brief  Test of Matrix class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_MatrixDense)

#include "TestDriver.hh"
#include "matrix/MatrixDense.hh"
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

//---------------------------------------------------------------------------//
// TEST DEFINITIONS
//---------------------------------------------------------------------------//

// Test of basic public interface
int test_MatrixDense(int argc, char *argv[])
{
  typedef MatrixDense Mat;
  typedef Vector Vec;

  // n * n
  {
    // Create test matrix
    int n = 5;
    Mat A(n, n, 0.0);
    TEST(A.number_rows()    == n);
    TEST(A.number_columns() == n);

    double row0[] = {-2.0, 1.1, 0.0, 0.0, 0.0};
    double row1[] = {1.0, -2.0, 1.1, 0.0, 0.0};
    double row2[] = {0.0, 1.0, -2.0, 1.1, 0.0};
    double row3[] = {0.0, 0.0, 1.0, -2.0, 1.1};
    double row4[] = {0.0, 0.0, 0.0, 1.0, -2.0};

    TEST(A.insert_row(0, row0, A.INSERT));
    TEST(A.insert_row(1, row1, A.INSERT));
    TEST(A.insert_row(2, row2, A.INSERT));
    TEST(A.insert_row(3, row3, A.INSERT));
    TEST(A.insert_row(4, row4, A.INSERT));
    A.assemble();
    A.print_matlab("A.out");

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
    Y.display();
    double ref[] =
    { 1.1, 2.4, 2.9, 3.6, -23.0 };
    for (int i = 0; i < Y.size(); i++)
    {
      Y.display();
      TEST(soft_equiv(Y[i], ref[i]));
    }

    // Transpose
    A.multiply_transpose(Y, X);
    double ref2[] =
    { 0.2, -0.69, 0.44, -27.01, 49.96 };
    for (int i = 0; i < X.size(); i++)
    {
      TEST(soft_equiv(X[i], ref2[i]));
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
    Mat A(m, n);
    TEST(A.number_rows()    == m);
    TEST(A.number_columns() == n);
    /*
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
      TEST(soft_equiv(Y[i], ref[i]));
    }

    // Transpose
    A.multiply_transpose(Y, X);

    double ref2[] =
    { 79.0, 100.0 };
    for (int i = 0; i < X.size(); i++)
    {
      TEST(soft_equiv(X[i], ref2[i]));
    }

    A.display();
    X.display();
    Y.display();

  } // end m * n

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_MatrixDense.cc
//---------------------------------------------------------------------------//
