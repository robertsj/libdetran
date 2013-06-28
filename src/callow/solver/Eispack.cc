//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Eispack.cc
 *  @brief Eispack member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Eispack.hh"
#include <complex>

namespace callow
{

//----------------------------------------------------------------------------//
Eispack::Eispack(const double tol,
                 const int    maxit)
  : EigenSolver(tol, maxit, "eispack")
{
  /* ... */
}


void Eispack::LZHES(MatrixDense &A, MatrixDense &B, MatrixDense &Z)
{
  Require(A.number_columns() == A.number_rows());
  Require(B.number_columns() == B.number_rows());
  Require(Z.number_columns() == Z.number_rows());
  Require(A.number_columns() == B.number_rows());
  Require(A.number_columns() == Z.number_rows());

  A.display();
  B.display();

  typedef std::complex<double> complex;
  using std::abs;
  using std::sqrt;

  int n = A.number_columns();
  if (n <= 1)
    return;

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      Z(i, j) = 0.0;

  DISP(" B to T ")

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
      double a = 0.0;
      double b = 0.0;
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
  if (n == 2)
    return;
  A.display();
  B.display();
  DISP(" A to hess ")
  // Reduce A to upper hessenberg, keeping B as triangle
  for (int k = 0; k < n - 2; ++k)
  {
    DISP(" k= " << k)
    int JM2 = n - 1 - k;
    int JP1 = k + 1;
    for (int ii = 0; ii < JM2; ++ii)
    {

      int I   = n - ii;
      int IM1 = I - 1;
      int IMJ = I - k;
      DISP(I << " " << IM1 <<" " << IMJ)
      double w = A(I, k);
      double z = A(IM1, k);
      double y = 0.0;
      // interchange rows if needed
      if (abs(w) >= abs(z))
      {
        DISP("interchange")

        for (int K = k; K < n; ++K)
        {
          y         = A(I, K);
          A(I, K)   = A(IM1, K);
          A(IM1, K) = y;
        }
        for (int K = IM1; K < n; ++K)
        {
          y         = B(I, K);
          B(I, K)   = B(IM1, K);
          B(IM1, K) = y;
        }
      }
      z = A(I, k);
      if (z != 0.0)
      {
        y = z / A(IM1, k);
        for (int K = JP1; K < n; ++K)
          A(I, K) = A(I, K) - y * A(IM1, K);
        for (int K = IM1; K < n; ++K)
          B(I, K) = B(I, K) - y * B(IM1, K);
      }
      w = B(I, IM1);
      z = B(I, I);
      if (abs(w) >= abs(z))
      {
        // MUST INTERCHANGE COLUMNS
        for (int K = 0; K < I; ++K)
        {
          y = B(K, I);
          B(K, I) = B(K, IM1);
          B(K, IM1) = y;
        }
        for (int K = 0; K < n; ++K)
        {
          y = A(K, I);
          A(K, I) = A(K, IM1);
          A(K, IM1) = y;
        }
        for (int K = IMJ; K < n; ++K)
        {
          y = Z(K, I);
          Z(K, I) = Z(K, IM1);
          Z(K, IM1) = y;
        }
      }
      z = B(I, IM1);
      if (z != 0.0)
      {
        y = z / B(I, I);
        for (int K = 0; K < IM1; ++K)
          B(K, IM1) = B(K, IM1) - y * B(K, I);
        B(I, IM1) = 0.0;
        for (int K = 0; K < n; ++K)
          A(K, IM1) = A(K, IM1) - y * A(K, I);

        for (int K = IMJ; K < n; ++K)
          Z(K, IM1) = Z(K, IM1) - y * Z(K, I);
      }
    }
    A(JP1 + 1, k) = 0.0;
  }

}

//----------------------------------------------------------------------------//
void Eispack::LZIT(MatrixDense &A,
                   MatrixDense &B,
                   MatrixDense &Z,
                   MatrixDense &V,
                   Vector      &L)
{




}

} // end namespace callow


//----------------------------------------------------------------------------//
//              end of file Eispack.cc
//----------------------------------------------------------------------------//


