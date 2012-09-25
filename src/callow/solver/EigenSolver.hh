//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   EigenSolver.hh
 *  @brief  EigenSolver
 *  @author Jeremy Roberts
 *  @date   Sep 24, 2012
 */
//---------------------------------------------------------------------------//

#ifndef EIGENSOLVER_HH_
#define EIGENSOLVER_HH_

#include "callow/callow_config.hh"
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
 *  We solve problems of the form
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
 *  Currently, we implement only the power method, inverse
 *  iteration, and Rayleigh quotient.  However, other solvers
 *  are available if SLEPc is enabled.  Additionally, our
 *  structure is really intended for the dominant mode.
 *  Other modes will require different handling.
 */
template <class T>
class EigenSolver
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  enum status
  {
    SUCCESS, MAXIT, DIVERGE
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef typename MatrixBase<T>::SP_matrix               SP_matrix;
  typedef typename LinearSolver<T>::SP_solver             SP_solver;
  typedef typename Vector<T>::SP_vector                   SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  EigenSolver(const double    tol = 1e-6,
              const int       maxit = 100,
              std::string     name = "solver")
    : d_tolerance(tol)
    , d_maximum_iterations(maxit)
    , d_residual_norm(maxit + 1, 0)
    , d_number_iterations(0)
    , d_monitor_level(1)
    , d_name(name)
  {
    Require(d_tolerance >= 0.0);
    Require(d_maximum_iterations > 0);
  }

  virtual ~EigenSolver(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Sets the operators for the linear system to solve.
   *
   *  This allows for the system
   *  @param A      right side operator
   *  @param B      left side operator
   */
  virtual void set_operators(SP_matrix A,
                             SP_matrix B = SP_matrix(0))
  {
    Require(A);
    d_A = A;
    Ensure(d_A->number_rows() == d_A->number_columns());
    if (B)
    {
      d_B = B;
      Ensure(d_B->number_rows() == d_B->number_columns());
      Ensure(d_B->number_rows() == d_A->number_columns());
      // create linear solver
    }
  }


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
   *  @param x    vector to fill with solution
   *  @param x0   initial guess
   */
  int solve(Vector<T> &x, Vector<T> &x0)
  {
    Require(x.size() == d_A->number_rows());
    if (x0.size())
      Require(x0.size() == d_A->number_rows());
    d_status = MAXIT;
    solve_impl(x, x0);
    if (d_status ==  MAXIT)
    {
      printf("*** %s did not converge within the maximum number of iterations\n",
             d_name.c_str());
    }
    return d_status;
  }

  /**
   *  @brief Solve the eigenvalue problem (SP variant)
   */
  int solve(SP_vector x, SP_vector x0)
  {
    return solve(*x, *x0);
  }

  /// return the residual norms
  std::vector<T> residual_norms()
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
  SP_solver linearsolver()
  {
    return d_solver;
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// convergence tolerance
  double d_tolerance;
  /// maximum iterations allowed
  int d_maximum_iterations;
  /// solver name
  std::string d_name;
  /// norm of successive residuals
  std::vector<T> d_residual_norm;
  /// number of iterations performed
  int d_number_iterations;
  /// left side operator
  SP_matrix d_A;
  /// right side operator (can be null)
  SP_matrix d_B;
  /// solver for inverting A
  SP_solver d_solver;
  /// diagnostic level
  int d_monitor_level;
  /// eigenvalue
  T d_lambda;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  // print out iteration and residual
  virtual bool monitor(int it, T l, T r)
  {
    // record the iteration and residual norm
    d_number_iterations = it;
    d_residual_norm[it] = r;
    // echo the residual
    if (d_monitor_level > 1)
      printf("iteration: %5i  eigenvalue: %12.8e    residual: %12.8e \n", it, l, r);
    // send a signal
    if (r < d_tolerance)
    {
      if (d_monitor_level > 0)
      {
        printf("*** %s converged in %5i iterations with a residual of %12.8e \n",
               d_name.c_str(), it, r);
      }
      d_status = SUCCESS;
      return true;
    }
    return false;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  virtual void solve_impl(Vector<T> &x, Vector<T> &x0) = 0;

private:

  int d_status;

};

} // end namespace callow

#endif // EIGENSOLVER_HH_ 

//---------------------------------------------------------------------------//
//              end of file Eigensolver.hh
//---------------------------------------------------------------------------//
