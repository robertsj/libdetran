//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenSolver.hh
 *  @brief EigenSolver class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_EIGENSOLVER_HH_
#define callow_EIGENSOLVER_HH_

#include "callow/callow_config.hh"
#include "callow/utils/CallowDefinitions.hh"
#include "LinearSolver.hh"
#include "callow/matrix/MatrixBase.hh"
#include "callow/vector/Vector.hh"
#include "utilities/SP.hh"
#include <cstdio>
#include <string>
#include <vector>

namespace callow
{

/**
 *  @class EigenSolver
 *  @brief Base class for iterative eigensolvers
 *
 *  We solve generalized eigenvalue problems of the form
 *  @f[
 *      \mathbf{A}x = \lambda \mathbf{B} x
 *  @f]
 *  via iterative methods.  For the case when
 *  @f$ \mathbf{B} \ne \mathbf{I} @f$ we need
 *  to solve
 *  @f[
 *      \mathbf{B}^{-1}\mathbf{A}x = \lambda x
 *  @f]
 *  A system is "solved" when the
 *  norm of the residual is small enough or some maximum
 *  iteration count is reached.  The residual is defined
 *  @f[
 *      \mathbf{A}x - \lambda \mathbf{B} x
 *  @f]
 *  and its norm is
 *  @f[
 *    r =  || (\mathbf{A} - \lambda \mathbf{B})x||
 *  @f]
 *  By default, the L2 norm is used, though L1 and Linf
 *  are also recorded can can be used. In some cases,
 *  a different norm is warranted, perhaps based on
 *  physics.  This can be implemented by derived classes.
 *
 *  Currently, we implement only the power method and
 *  nonlinear Arnoldi method.  Additionally, our
 *  structure is really intended for the dominant mode.
 *  Other modes will require different handling.
 */
class CALLOW_EXPORT EigenSolver
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<EigenSolver>       SP_solver;
  typedef MatrixBase::SP_matrix                   SP_matrix;
  typedef LinearSolver::SP_solver                 SP_linearsolver;
  typedef Preconditioner::SP_preconditioner       SP_preconditioner;
  typedef Vector::SP_vector                       SP_vector;
  typedef detran_utilities::InputDB::SP_input     SP_db;
  typedef detran_utilities::size_t                size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *   @brief Constructor
   *   @param     tol     tolerance
   *   @param     maxit   maximum number of iterations
   *   @param     name    solver name
   */
  EigenSolver(const double    tol   = 1e-6,
              const int       maxit = 100,
              std::string     name  = "solver");

  /// Destructor
  virtual ~EigenSolver(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Sets the operators for the problem.
   *
   *  This allows for the system
   *  @param A      left side operator
   *  @param B      optional right side operator (to be inverted)
   *  @param db     optional database for solver and preconditioner options
   */
  virtual void set_operators(SP_matrix  A,
                             SP_matrix  B  = SP_matrix(0),
                             SP_db      db = SP_db(0));


  /// Set the preconditioner for a generalized eigenvalue problem
  virtual void set_preconditioner(SP_preconditioner P,
                                  const int         side = LinearSolver::LEFT);

  /// Set the preconditioner matrix for a generalized eigenvalue problem
  virtual void set_preconditioner_matrix(SP_matrix P,
                                         const int side = LinearSolver::LEFT);

  /**
   *  @brief Set the convergence criteria
   *  @param tol  m tolerance (||r_n|| < tol)
   *  @param maxit  maximum iterations (n < maxit)
   */
  void set_tolerances(const double tol, const int maxit)
  {
    d_tolerance = tol;
    d_maximum_iterations = maxit;
    Require(d_tolerance > 0.0);
    Require(d_maximum_iterations >= 0);
  }

  /**
   *  @brief Set the diagnostic level.
   *
   *  Higher levels means more output.
   *
   *  @param v  diagnostic level
   */
  void set_monitor_level(const int v)
  {
    d_monitor_level = v;
  }

  /**
   *  @brief Solve the eigenvalue problem
   *
   *  Upon return, the initial guess is *not*
   *  guaranteed to be unchanged, as it may be used
   *  as temporary storage.
   *
   *  @param x    vector to fill with solution
   *  @param x0   initial guess
   */
  int solve(Vector &x, Vector &x0);

  /**
   *  @brief Solve the eigenvalue problem (SP variant)
   */
  int solve(SP_vector x, SP_vector x0)
  {
    return solve(*x, *x0);
  }

  /// Return the residual norms
  std::vector<double> residual_norms()
  {
    return d_residual_norm;
  }

  /// Get the left operator
  SP_matrix A()
  {
    return d_A;
  }

  /// Get the right operator
  SP_matrix B()
  {
    return d_B;
  }

  /// Get the linear solver
  SP_linearsolver linearsolver()
  {
    return d_solver;
  }

  /// Get the eigenvalue
  double eigenvalue()
  {
    return d_lambda;
  }

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// convergence tolerance
  double d_tolerance;
  /// maximum iterations allowed
  int d_maximum_iterations;
  /// solver name
  std::string d_name;
  /// norm of successive residuals
  std::vector<double> d_residual_norm;
  /// number of iterations performed
  int d_number_iterations;
  /// left side operator
  SP_matrix d_A;
  /// right side operator (can be null)
  SP_matrix d_B;
  /// solver for inverting A
  SP_linearsolver d_solver;
  /// diagnostic level
  int d_monitor_level;
  /// eigenvalue
  double d_lambda;
  /// solver status
  int d_status;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /**
   *  @brief Print out iteration and residual
   *
   *  Monitoring can be overloaded so that method-specific criteria
   *  can be employed.
   *
   *  @param it iteration
   *  @param l  eigenvalue estimate (dominant, if others computed)
   *  @param r  residual
   */
  virtual bool monitor(int it, double l, double r);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  virtual void solve_impl(Vector &x, Vector &x0) = 0;

};

CALLOW_TEMPLATE_EXPORT(detran_utilities::SP<EigenSolver>)

} // end namespace callow

#include "EigenSolver.i.hh"

#endif // callow_EIGENSOLVER_HH_

//----------------------------------------------------------------------------//
//              end of file Eigensolver.hh
//----------------------------------------------------------------------------//
