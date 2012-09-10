//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_OperatorMatrix.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of OperatorMatrix class.
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                 \
        FUNC(test_OperatorMatrix)

// Detran test
#include "utilities/TestDriver.hh"
#include "solvers/OperatorMatrix.hh"
#include <iostream>

// Setup
#include "solvers/test/operator_fixture.hh"

using namespace detran_test;
using namespace detran;
using namespace detran_utilities;
using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_OperatorMatrix_actual();

// Test of basic public interface
int test_OperatorMatrix(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Call actual test.
  int result = test_OperatorMatrix_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;
}

// Test of basic public interface
int test_OperatorMatrix_actual()
{
  typedef detran_utilities::size_t size_t;

  // Create test matrix
  size_t  n = 5;
  TestOperatorMatrix A(n);
  TEST(A.number_rows() == n);
  TEST(A.number_columns() == n);

  // Create two vectors
  Vector X(n);
  Vector Y(n);

  for (int i = 0; i < X.size(); i++)
  {
    double value = i*i;
    X.insert_values(1, &i, &value);
  }
  X.assemble();
  Y.assemble();

  // Multiply
  A.multiply(X, Y);

  // Test the vector output.
  double ref[] = {1, 2, 2, 2, -23};
  for (int i = 0; i < X.size(); i++)
  {
    TEST(detran_utilities::soft_equiv(Y[i], ref[i]));
  }

  A.display();
  X.display();
  Y.display();

  Y[0] = 2.0;
  Y.display();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Operator.cc
//---------------------------------------------------------------------------//
