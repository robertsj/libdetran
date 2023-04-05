//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Eispack.cc
 *  @brief Test of Eispack class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "callow/utils/Initialization.hh"
#include "callow/solver/Eispack.hh"
#include <iostream>

using namespace callow;

using std::cout;
using std::endl;

#define DISP(c) cout << c << endl;
#define W() std::cout << "line " << __LINE__ << " on file " << __FILE__ << std::endl;

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int n = 5;
double a[] =
{ 0.6958, 0.2468, 0.2228, 0.0258, 0.6501, 0.9880, 0.9525, 0.0023, 0.6045,
  0.2479, 0.9671, 0.4849, 0.4617, 0.3953, 0.5141, 0.3681, 0.3422, 0.3677,
  0.4179, 0.6940, 0.3904, 0.1448, 0.9227, 0.1692, 0.4964 };

MatrixDense::SP_matrix get_A()
{
  MatrixDense::SP_matrix A(new MatrixDense(n, n, 0.0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      (*A)(i, j) = a[j + i * n] + i;
  A->print_matlab("A.out");
  return A;
}

MatrixDense::SP_matrix get_B()
{
  MatrixDense::SP_matrix B(new MatrixDense(n, n, 0.0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      (*B)(i, j) = a[j + i * n];
  B->print_matlab("B.out");
  return B;
}

//----------------------------------------------------------------------------//
TEST(Eispack, Standard)
{
  Eispack S(1.e-5, 1000);
  auto A = get_A();
  Require(A);
  S.set_operators(A);
  // test single eigenvalue and vector
  {
    Vector X0(n, 0.0), X(n, 0.0);
    S.solve(X, X0);
    EXPECT_NEAR(S.eigenvalue(), 1.216310176599158e+01, 1.0e-12);
    double ref[] = {-5.332899955518561e-02, -2.210866243722326e-01,
        -3.907759112294125e-01, -5.515953428318574e-01, -7.009375773199511e-01};
    for (int i = 0; i < n; ++i)
    {
      EXPECT_NEAR(X[i], ref[i], 1.0e-12);
    }
  }

  // test all eigenvalues
  {
    MatrixDense V_R(n, n, 0.0), V_I(n, n, 0.0);
    Vector      E_R(n, 0.0),    E_I(n, 0.0);
    S.solve_complete(V_R, V_I, E_R, E_I);
    double ref_r[] = {1.216310176599158e+01, -5.783293144882002e-02,
        -5.783293144882002e-02, 4.665470070517458e-01, 5.103170898543222e-01};
    double ref_i[] = {0, 3.504177633722163e-01, -3.504177633722163e-01, 0, 0};
    for (int i = 0; i < n; ++i)
    {
      EXPECT_NEAR(E_R[i], ref_r[i], 1.0e-12);
      EXPECT_NEAR(E_I[i], ref_i[i], 1.0e-12);
    }
  }
}

//----------------------------------------------------------------------------//
TEST(Eispack, General)
{
  Eispack S(1.e-5, 1000);
  S.set_operators(get_A(), get_B());

  // test single eigenvalue and vector
  {
    Vector X0(n, 0.0), X(n, 0.0);
    S.solve(X, X0);
    EXPECT_NEAR(S.eigenvalue(), 8.151729277070530e+00, 1.0e-12);
    double ref[] = {3.324109144209688e-01, -6.153722684382446e-02,
        -6.525478413205275e-01, -6.731740329299938e-01, -8.206210978810474e-02};
    X.display();
    for (int i = 0; i < n; ++i)
    {
       EXPECT_NEAR(X[i], ref[i], 1.0e-12);
    }
  }

  // test all eigenvalues
  {
    MatrixDense V_R(n, n, 0.0), V_I(n, n, 0.0);
    Vector      E_R(n, 0.0),    E_I(n, 0.0);
    S.solve_complete(V_R, V_I, E_R, E_I); E_R.display(); E_I.display();
    double ref_r[] = {8.151729277070530e+00, 1.0, 1.0, 1.0, 1.0};
    double ref_i[] = {0, 0, 0, 0, 0};
    E_R.display();
    for (int i = 0; i < n; ++i)
    {
      EXPECT_NEAR(E_R[i], ref_r[i], 1.0e-12);
      EXPECT_NEAR(E_I[i], ref_i[i], 1.0e-12);
    }
  }
}

//----------------------------------------------------------------------------//
//              end of test_Eispack.cc
//----------------------------------------------------------------------------//
