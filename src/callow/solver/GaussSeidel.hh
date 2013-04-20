//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussSeidel.hh
 *  @author robertsj
 *  @date   Sep 14, 2012
 *  @brief  GaussSeidel class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_GAUSSSEIDEL_HH_
#define callow_GAUSSSEIDEL_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  @class GaussSeidel
 *  @brief Uses Gauss-Seidel iteration to solve a system
 *
 *  Gauss-Seidel iteration in matrix form uses a splitting of
 *  the form
 *  @f[
 *      \mathbf{A} = \mathbf{L} + \mathbf{U} +  \mathbf{D} \, ,
 *  @f]
 *  which are strictly lower and upper triangle and diagonal,
 *  respectively.  The Gauss-Seidel iteration is then
 *  @f[
 *      \mathbf{D + L} x^{n+1} = -\mathbf{U}x^{n} + b
 *  @f]
 *  or
 *  @f[
 *      x^{n+1} = \overbrace{-\mathbf{D+L}^{-1}(\mathbf{U})}^
 *                          {\mathbf{M}}x^{n} + \mathbf{D+L}^{-1}b \, .
 *  @f]
 *  Alternatively, one can swap @f$ U @f$ and @f$ L @f$ to produce
 *  the backward Gauss-Seidel iteration.  If used together, one
 *  has symmetric Gauss-Seidel iteration.
 *
 *  It can be shown that Gauss-Seidel converges if \ref Jacobi
 *  converges.
 *
 *  Because we use a sparse matrix, we actually access the elements
 *  directly rather than via indexing.
 *
 */
class GaussSeidel: public LinearSolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  GaussSeidel(const double atol,
              const double rtol,
              const int maxit,
              const double omega = 1.0,
              bool successive_norm = false);

  virtual ~GaussSeidel(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver::d_absolute_tolerance;
  using LinearSolver::d_relative_tolerance;
  using LinearSolver::d_maximum_iterations;
  using LinearSolver::d_residual;
  using LinearSolver::d_number_iterations;
  using LinearSolver::d_A;
  using LinearSolver::d_P;
  using LinearSolver::d_norm_type;

  // Relaxation parameter
  double d_omega;

  // Use norm of residual (false=default) or successive iterates
  bool d_successive_norm;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  @param b  right hand side
   *  @param x  unknown vector
   */
  void solve_impl(const Vector &b, Vector &x);

};

} // end namespace callow

// Inline member definitions
#include "GaussSeidel.i.hh"


#endif /* callow_GAUSSSEIDEL_HH_ */
