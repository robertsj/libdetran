//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Eispack.i.hh
 *  @brief Eispack inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_EISPACK_I_HH_
#define callow_EISPACK_I_HH_

#define DISP(c) std::cout << c << std::endl;

namespace callow
{

//----------------------------------------------------------------------------//
inline void Eispack::solve_impl(Vector &x, Vector &x0)
{

}

//----------------------------------------------------------------------------//
inline void Eispack::qzhes(MatrixDense &A, MatrixDense &B, MatrixDense &Z)
{
  Require(A.number_columns() == A.number_rows());
  Require(B.number_columns() == B.number_rows());
  Require(Z.number_columns() == Z.number_rows());
  Require(A.number_columns() == B.number_rows());
  Require(A.number_columns() == Z.number_rows());

  using std::abs;
  using std::sqrt;

  int n = A.number_columns();
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      Z(i, j) = 0.0;

  if (n == 1) return;

  //------------------------------------------------------------------------//
  // TRANSFORM GENERAL B TO UPPER TRIANGLE T VIA HOUSEHOLDER REFLECTIONS
  //------------------------------------------------------------------------//
  DISP("Transforming B to T...")
  for (int k = 0; k < n - 1; ++k)
  {
    DISP(" column=" << k)
    /*
     *  Let X[0:k] = B[0:k, k].  We zero out X[0:k-1] via Householder:
     *    X <-- P*X = (I - 2*VV'/V'V)*X ,
     *  where
     *    V = -sign(X(k))*||X||*E_k - X ,
     *  E_k has all zeros except a one in the kth location, and ||.|| denotes
     *  the L2 norm.
     */
    double X_L2  = 0.0;
    for (int i = k + 1; i < n; ++i)
      X_L2 += pow(B(i, k), 2);
    DISP(" normalizing X " << k << B(k, k))

    if (X_L2 == 0.0) break;
    B(k, k) += sign(sqrt(X_L2 + pow(B(k, k), 2)), B(k, k));
    DISP(" normalizing X " << k)

    X_L2 = sqrt(X_L2+pow(B(k, k), 2));
    for (int i = k; i < n; ++i)
      B(i, k) /= X_L2;
    DISP(" A=" << k)

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

  //------------------------------------------------------------------------//
  // REDUCE A TO UPPER HESSENBERG, KEEPING B UPPER TRIANGLE
  //------------------------------------------------------------------------//

  /*
   *  We have to do very nearly the same thing to reduce A.  For each
   *  partial column we eliminate, we must apply a second orthogonal
   *  operator to restore B's triangular structure.
   *
   */

  if (n == 2) return;

  for (int k = 0; k < n - 2; ++k)
  {
    DISP(" below column=" << k)
    //     .......... for l=n-1 step -1 until k+1 do -- ..........
    for (int lb = 0; lb < n - k - 1; ++lb)
    {
      int l  = n - lb;
      int l1 = l + 1;
      DISP("lb=" << lb << " l=" << l << " l1=" << l1)
      //    .......... zero A(l+1, k) ..........
      double s = std::abs(A(l, k)) + std::abs(A(l1, k));
      if (s == 0.0) break;
      double u1 = A(l, k) / s;
      double u2 = A(l1, k) / s;
      double  r = sign(std::sqrt(u1 * u1 + u2 * u2), u1);
      double v1 = -(u1 + r) / r;
      double v2 = -u2 / r;
      u2        = v2 / v1;
      double  t = 0.0;

      for (int j = k; j < n; ++j)
      {
        t = A(l, j) + u2 * A(l1, j);
        A(l, j) = A(l, j) + t * v1;
        A(l1, j) = A(l1, j) + t * v2;
      }

      A(l1, k) = 0.0;

      for (int j = l; j < n; ++j)
      {
        t = B(l, j) + u2 * B(l1, j);
        B(l, j)  = B(l, j) + t * v1;
        B(l1, j) = B(l1, j) + t * v2;
      }
      //     .......... zero B(l+1,l) ..........
      s = std::abs(B(l1, l1)) + std::abs(B(l1, l));
      if (s == 0.0) break;
      u1 = B(l1, l1) / s;
      u2 = B(l1, l) / s;
      r = sign(std::sqrt(u1 * u1 + u2 * u2), u1);
      v1 = -(u1 + r) / r;
      v2 = -u2 / r;
      u2 = v2 / v1;

      for (int i = 0; i < l1; ++i)
      {
        t         = B(i, l1) + u2 * B(i, l);
        B(i, l1) += t * v1;
        B(i, l)  += t * v2;
      }

      B(l1, l) = 0.0;

      for (int i = 0; i < n; ++i)
      {
        t         = A(i, l1) + u2 * A(i, l);
        A(i, l1) += t * v1;
        A(i, l)  += t * v2;
      }

      for (int i = 0; i < n; ++i)
      {
        t         = Z(i, l1) + u2 * Z(i, l);
        Z(i, l1) += t * v1;
        Z(i, l)  += t * v2;
      }
    }
  }
}

} // end namespace callow

#endif /* callow_EISPACK_I_HH_ */

//----------------------------------------------------------------------------//
//              end of file Eispack.i.hh
//----------------------------------------------------------------------------//
