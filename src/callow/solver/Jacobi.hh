//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Jacobi.hh
 *  @brief Jacobi class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_JACOBI_HH_
#define callow_JACOBI_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  @class Jacobi
 *  @brief Uses (weighted) Jacobi iteration to solve a system
 *
 *  Jacobi iteration in matrix form uses a splitting of
 *  the form
 *  @f[
 *      \mathbf{A} = \mathbf{L} + \mathbf{U} +  \mathbf{D} \, ,
 *  @f]
 *  which are strictly lower and upper triangle and diagonal,
 *  respectively.  The Jacobi iteration is then
 *  @f[
 *      \mathbf{D} x^{n+1} = -\omega (\mathbf{L}+\mathbf{U})x^{n} +
 *                           (1-\omega) b
 *  @f]
 *  or
 *  @f[
 *      x^{n+1} = \overbrace{-\omega\mathbf{D}^{-1}(\mathbf{L}+\mathbf{U})}^
 *                          {\mathbf{M}}x^{n} + \omega\mathbf{D}^{-1}b \, .
 *  @f]
 *  The procedure converges if the iteration matrix is bounded, i.e.
 *  @f[
 *     \rho(\mathbf{M}) < 1 \, ,
 *  @f]
 *  which (for @f$ \omega = 1 @f$) is guaranteed if \f$ \mathbf{A} \f$ is
 *  strictly diagonally dominant, meaning that
 *  @f[
 *      |a_{ii}| > \sum_{j \neq i} |a_{ij}| \, .
 *  @f]
 */
class Jacobi: public LinearSolver
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef LinearSolver Base;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  Jacobi(const double atol,
         const double rtol,
         const int    maxit,
         const double omega = 1.0,
         bool         successive_norm = false);

  virtual ~Jacobi(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Use norm of residual (false=default) or successive iterates
  bool d_successive_norm;
  /// Weighting factor
  double d_omega;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /**
   *  @param b  right hand side
   *  @param x  unknown vector
   */
  void solve_impl(const Vector &b, Vector &x);

};

} // end namespace callow

// Inline member definitions
#include "Jacobi.i.hh"

#endif // callow_JACOBI_HH_

//----------------------------------------------------------------------------//
//              end of file Jacobi.hh
//----------------------------------------------------------------------------//
