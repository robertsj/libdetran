//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Eispack.cc
 *  @brief Test of Eispack class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//
// LIST OF TEST FUNCTIONS
#define TEST_LIST                 \
        FUNC(test_Eispack)

#include "utilities/TestDriver.hh"
#include "callow/utils/Initialization.hh"
#include "callow/solver/Eispack.hh"
#include <iostream>

using namespace callow;
using namespace detran_test;
using detran_utilities::soft_equiv;
using std::cout;
using std::endl;

#define DISP(c) cout << c << endl;

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
int test_Eispack(int argc, char *argv[])
{
  int n = 5;
  double a[] =
  { 0.6958, 0.2468, 0.2228, 0.0258, 0.6501, 0.9880, 0.9525, 0.0023, 0.6045,
    0.2479, 0.9671, 0.4849, 0.4617, 0.3953, 0.5141, 0.3681, 0.3422, 0.3677,
    0.4179, 0.6940, 0.3904, 0.1448, 0.9227, 0.1692, 0.4964 };
  MatrixDense::SP_matrix A(new MatrixDense(n, n, 0.0));
  MatrixDense::SP_matrix B(new MatrixDense(n, n, 0.0));
  MatrixDense::SP_matrix Z(new MatrixDense(n, n, 0.0));

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      (*A)(i, j) = a[j + i * n] + i;
      (*B)(i, j) = a[j + i * n];
    }
  }
  //MatrixDense C(A), D(B), X(n, n, 0.0);
  DISP("B")
  B->print_matlab("B.out");
  A->print_matlab("A.out");

  Eispack S(1.e-5, 1000);
  S.set_operators(A);

  MatrixDense V_R(n, n, 0.0), V_I(n, n, 0.0);
  Vector E_R(n, 0.0), E_I(n, 0.0);

  S.solve_complete(V_R, V_I, E_R, E_I);
  E_R.display();
  E_I.display();

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Eispack.cc
//----------------------------------------------------------------------------//
