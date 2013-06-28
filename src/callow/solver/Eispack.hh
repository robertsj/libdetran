//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Eispack.hh
 *  @brief Eispack class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_EISPACK_HH_
#define callow_EISPACK_HH_

#include "EigenSolver.hh"
#include "matrix/MatrixDense.hh"
#include <cmath>
#include <limits>

namespace callow
{

/**
 *  @class Eispack
 *  @brief Solve a dense eigenvalue problem with the QZ method
 *
 *  The basic algorithm was influence from the old Eispack software.  Some
 *  things have been improved and adapted as needed.
 */

class Eispack: public EigenSolver
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  Eispack(const double tol   = 1e-6,
          const int    maxit = 100);

  virtual ~Eispack(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Dense matrices A and B to upper hessenberg H and upper triangle T
   *
   *  Householder reflections are used to triangulate B.
   *
   *  @param A    right hand matrix to become upper hessenberg
   *  @param B    left hand matrix to become upper triangle
   *  @param Z    contains the product of the right hand transformations, needed
   *              for producing the vectors
   */
  static void qzhes(MatrixDense &A, MatrixDense &B, MatrixDense &Z);


  static void LZHES(MatrixDense &A, MatrixDense &B, MatrixDense &X);
  static void  LZIT(MatrixDense &A,
                    MatrixDense &B,
                    MatrixDense &X,
                    MatrixDense &V,
                    Vector      &L);
private:

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  void solve_impl(Vector &x, Vector &x0);

};

/// Return value of a with sign of b
inline double sign(const double a, const double b)
{
  double s = (b < 0) ? -1 : 1;
  return std::abs(a) * s;
}

/// Estimate unit roundoff in quantities of size x
inline double epslon(const double x)
{
  return std::numeric_limits<double>::epsilon() * std::abs(x);
}


} // end namespace callow

#include "Eispack.i.hh"

#endif /* callow_EISPACK_HH_ */

//----------------------------------------------------------------------------//
//              end of file Eispack.hh
//----------------------------------------------------------------------------//
