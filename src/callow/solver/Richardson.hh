//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Richardson.hh
 *  @author robertsj
 *  @date   Sep 13, 2012
 *  @brief  Richardson class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_RICHARDSON_HH_
#define callow_RICHARDSON_HH_

#include "LinearSolver.hh"

namespace callow
{

/**
 *  @class Richardson
 *  @brief Uses (modified) Richardson iteration to solve a system
 *
 *  Richardson iteration solves a linear system via the
 *  process
 *  @f[
 *     x^{(n+1)} = (\mathbf{I - \omega A})x^{(n)} + \omega b
 *  @f]
 *  where \f$ \omega \f$ is something like a relaxation factor
 *  that takes on values between (roughly) 0 and 2.  By default,
 *  \f$ \omega = 1 \f$.
 *
 */
class Richardson: public LinearSolver
{

public:

  typedef LinearSolver Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  Richardson(const double atol, const double rtol, const int maxit,
             const double omega = 1.0);

  virtual ~Richardson(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  void set_omega(const double omega)
  {
    d_omega = omega;
  }

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

  /// relaxation factor
  double d_omega;

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
#include "Richardson.i.hh"

#endif /* callow_RICHARDSON_HH_ */
