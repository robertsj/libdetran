//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_OperatorShell.cc
 * \author Jeremy Roberts
 * \date   Aug 19, 2012
 * \brief  Test of OperatorShell class.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST         \
        FUNC(test_OperatorShell)

// Detran test
#include "TestDriver.hh"

#include "OperatorShell.hh"

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

int test_OperatorShell_actual();

// Test of basic public interface
int test_OperatorShell(int argc, char *argv[])
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Call actual test.
  int result = test_OperatorShell_actual();

  // Finalize PETSc
  PetscFinalize();

  return result;
}

// Test of basic public interface
int test_OperatorShell_actual()
{

  // Create OperatorShell
  detran_utilities::size_t n = 5;
  TestOperatorShell A(n);
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
  double ref[] = {1, 2, -X.size() * 5 + 2};
  for (int i = 0; i < X.size(); i++)
  {
    double ref = 2.0;
    if (i == 0) ref = 1.0;
    if (i == n - 1) ref = 2 - 25 * X.size() * X.size();
    TEST(detran_utilities::soft_equiv(Y[i], ref));
  }

  A.display();
  X.display();
  Y.display();

  Y[0] = 2.0;
  Y.display();
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_OperatorShell.cc
//---------------------------------------------------------------------------//
