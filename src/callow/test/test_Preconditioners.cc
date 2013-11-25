//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Preconditioners.cc
 *  @brief Test of Matrix class
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST            \
        FUNC(test_PCJacobi)  \
        FUNC(test_PCILU0)

#include "TestDriver.hh"
#include "preconditioner/PCJacobi.hh"
#include "preconditioner/PCILU0.hh"

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
//              end of test_Preconditioners.cc
//----------------------------------------------------------------------------//
