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
//#define DISP(c) cout << c << endl;

int main(int argc, char *argv[])
{
  callow_initialize(argc, argv);
  RUN(argc, argv);
  callow_finalize();
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

inline void triQ(MatrixDense &A, MatrixDense &B)
{
  int n = A.number_columns();
  // Reduce B to upper triangle
  for (int k = 0; k < n - 1; ++k)
  {
    double X_L2 = 0.0;
    for (int i = k + 1; i < n; ++i)
      X_L2 += pow(B(i, k), 2);
    if (X_L2 == 0.0)
      continue;
    B(k, k) += sign(sqrt(X_L2 + pow(B(k, k), 2)), B(k, k));
    X_L2 = sqrt(X_L2+pow(B(k, k), 2));
    for (int i = k; i < n; ++i)
      B(i, k) /= X_L2;
    for (int j = k; j < n; ++j)
    {
      double a = 0.0, b = 0.0;
      for (int i = k; i < n; ++i)
      {
        a += B(i, k) * A(i, j);
        b += B(i, k) * B(i, j);
      }
      a *= -2.0;
      b *= -2.0;
      for (int i = k; i < n; ++i)
      {
        A(i, j) += a * B(i, k);
        B(i, j) += b * B(i, k);
      }
    }
  }
}

/// Reduce B to T
inline void tri(MatrixDense &A, MatrixDense &B)
{
  int N = A.number_rows();
  int NM1 = N - 1;
  int II = 0;
  // loop over all but last column
  for (int I = 0; I < NM1; ++I)
  {
    double D = 0.0;
    double Y = 0.0;
    int IP1  = I + 1;
    // pick out the maximum element in this column (below the diagonal)
    for (int K = IP1; K < N; ++K)
    {
      Y = std::abs(B(K, I));
      if (Y < D) continue;
      D  = Y;
      II = K;
    }
    // jump to next column if this one already has zeros under the diagonal
    if (D == 0.0) continue;
    Y = B(I, I);
    // exchange rows if a subdiagonal is larger
    if (D >= Y)
    {
      for (int J = 0; J < N; ++J)
      {
        Y        = A(I,  J);
        A(I,  J) = A(II, J);
        A(II, J) = Y;
      }
      for (int J = I; J < N; ++J)
      {
        Y        = B(I,  J);
        B(I,  J) = B(II, J);
        B(II, J) = Y;
      }
    }
    for (int J = IP1; J < N; ++J)
    {
      double Y = B(J, I) / B(I, I);
      if (Y == 0.0) continue;
      for (int K = 0; K < N; ++K)
        A(J, K) = A(J, K) - Y * A(I, K);
      for (int K = IP1; K < N; ++K)
        B(J, K) = B(J, K) - Y * B(I, K);
      B(J, I) = 0.0;
    }
  }
}

// REDUCE A TO UPPER HESSENBERG FORM
inline void hess(MatrixDense &A, MatrixDense &B, MatrixDense &X)
{
  int N = A.number_rows();
  int NM1 = N - 1;
  int NM2 = N - 2; // 110
  if (NM2 < 1)
    return;
  int II = 0;
  double Y = 0;

  // INITIALIZE X
  for (int I = 0; I < N; ++I)
  {
    for (int J = 0; J < N; ++J)
    {
      X(I, J) = 0.0;
    }
    X(I, I) = 1.0;
  }

  for (int J = 0; J < NM2; ++J)
  {
    //DISP("J=" << J)
    int JM2 = NM1 - J;
    int JP1 = J + 1;
    for (II = 0; II < JM2; ++II)
    {
      //DISP("II=" << II)
      int I = N - II - 1;
      int IM1 = I - 1;
      int IMJ = I - J;
//      DISP("I=" << I)
//      DISP("IM1=" << IM1)
//      DISP("IMJ" << IMJ)
      double W = A(I, J);
      double Z = A(IM1, J);
      // MUST INTERCHANGE ROWS
      if (std::abs(W) >= std::abs(Z))
      {
        //DISP("interchange")
        for (int K = J; K < N; ++K)
        {
          Y = A(I, K);
          A(I, K) = A(IM1, K);
          A(IM1, K) = Y;
        } // 120
        for (int K = IM1; K < N; ++K)
        {
          Y = B(I, K);
          B(I, K) = B(IM1, K);
          B(IM1, K) = Y;
        } // 130
      }
      DISP("lala")
      Z = A(I, J); // 140
      if (Z != 0.0)
      {
        Y = Z / A(IM1, J);
        for (int K = JP1; K < N; ++K)
        {
          A(I, K) = A(I, K) - Y * A(IM1, K);
        } // 150
        for (int K = IM1; K < N; ++K)
        {
          B(I, K) = B(I, K) - Y * B(IM1, K);
        } // 160
      }
      // TRANSFORMATION FROM THE RIGHT
      W = B(I, IM1); // 170
      Z = B(I, I);
      if (std::abs(W) > std::abs(Z))
      {
        // MUST INTERCHANGE COLUMNS
        for (int K = 0; K < I; ++K)
        {
          Y = B(K, I);
          B(K, I) = B(K, IM1);
          B(K, IM1) = Y;
        } // 180
        for (int K = 0; K < N; ++K)
        {
          Y = A(K, I);
          A(K, I) = A(K, IM1);
          A(K, IM1) = Y;
        } // 190
        for (int K = IMJ; K < N; ++K)
        {
          Y = X(K, I);
          X(K, I) = X(K, IM1);
          X(K, IM1) = Y;
        }
      }
      Z = B(I, IM1); // 210
      if (Z != 0.0)
      {
        Y = Z / B(I, I);
        for (int K = 0; K < IM1; ++K)
        {
          B(K, IM1) = B(K, IM1) - Y * B(K, I);
        } // 220
        B(I, IM1) = 0.0;
        for (int K = 0; K < N; ++K)
        {
          A(K, IM1) = A(K, IM1) - Y * A(K, I);
        } // 230
        for (int K = IMJ; K < N; ++K)
        {
          X(K, IM1) = X(K, IM1) - Y * X(K, I);
        } // 240
      }
    } // 250
    A(JP1 + 1, J) = 0.0;
  }
  return; // 270
}

//----------------------------------------------------------------------------//
int test_Eispack(int argc, char *argv[])
{
  int n = 5;
  double a[] =
  { 0.6958, 0.2468, 0.2228, 0.0258, 0.6501, 0.9880, 0.9525, 0.0023, 0.6045,
    0.2479, 0.9671, 0.4849, 0.4617, 0.3953, 0.5141, 0.3681, 0.3422, 0.3677,
    0.4179, 0.6940, 0.3904, 0.1448, 0.9227, 0.1692, 0.4964 };
  MatrixDense A(n, n, 0.0), B(n, n, 0.0), Z(n, n, 0);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      A(i, j) = a[j + i * n] + i;
      B(i, j) = a[j + i * n];
    }
  }
  MatrixDense C(A), D(B), X(n, n, 0.0);
  DISP("B")
  B.display();
  B.print_matlab("B.out"); A.print_matlab("A.out");
  tri(A, B);
  DISP("A")
  A.display();
  DISP("B")
  B.display();
  DISP("");
  hess(A, B, X);
  A.display();
  B.display();
//  triQ(C, D);
//  DISP("C")
//  C.display();
//  DISP("D")
//  D.display();
//
//  Eispack S(1.e-5, 1000);
//  Eispack::LZHES(A, B, Z);

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_Eispack.cc
//----------------------------------------------------------------------------//
