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

namespace callow
{

/**
 *  @class Eispack
 *  @brief Solve a dense eigenvalue problem with either QR or QZ
 *
 *  This class is just a wrapper around old Eispack functionality.  We could
 *  link to LAPACK, but it's nice to have it all self-contained.
 */

class Eispack: public EigenSolver
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  Eispack(const double tol   = 1e-6,
          const int    maxit = 100,
          int          which = 1);

  virtual ~Eispack(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Sets the operators for the problem.  Eispack requires dense matrices.
  void set_operators(SP_matrix  A,
                     SP_matrix  B  = SP_matrix(0),
                     SP_db      db = SP_db(0));

  /// Solve for the complete eigenspectrum
  void solve_complete(MatrixDense &V_R,
                      MatrixDense &V_I,
                      Vector      &E_R,
                      Vector      &E_I);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Indicates largest, smallest, etc.
  int d_which_value;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  void solve_impl(Vector &x, Vector &x0);

};

extern "C"
{
// SUBROUTINE rg (nm, n, a, wr, wi, matz, z, iv1, fv1, ierr)
void rg_(const int* nm, const int* n, double* a, double* wr, double* wi,
         const int* matz, double* z, int* iv1, double* fv1, int* ierr);

// SUBROUTINE rgg (nm, n, a, b, alfr, alfi, beta, matz, z, ierr)
void rgg_(int*, int*, double*, double*, double*, double*, double*, int*,
          double*, int*);

// SUBROUTINE qzhes (nm, n, a, b, matz, z)
void qzhes_(int*, int*, double*, double*, bool*, double*);

// SUBROUTINE qzit (nm, n, a, b, eps1, matz, z, ierr)
void qzit_(int*, int*, double*, double*, double*, bool*, double*, int*);

// SUBROUTINE qzval (nm, n, a, b, alfr, alfi, beta, matz, z)
void qzval_(const int* nm, const int* n, double* a, double* b, double* alfr,
            double* alfi, double* beta, const bool* matz, double* z);



}

} // end namespace callow

#endif /* callow_EISPACK_HH_ */

//----------------------------------------------------------------------------//
//              end of file Eispack.hh
//----------------------------------------------------------------------------//
