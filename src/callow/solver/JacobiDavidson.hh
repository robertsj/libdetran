//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  JacobiDavidson.hh
 *  @brief JacobiDavidson class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef callow_JACOBIDAVIDSON_HH_
#define callow_JACOBIDAVIDSON_HH_

#include "EigenSolver.hh"
#include "callow/preconditioner/PCShell.hh"
#include "callow/matrix/MatrixShell.hh"

namespace callow
{

/**
 *  @class JacobiDavidson
 *  @brief Solve the eigenvalue problem with Jacobi-Davidson iteration.
 *
 *  The basic idea of J-D iteration is build a subspace based on
 *  an approximation to the original matrix (part of the "preconditioner").
 *  If that approximation is cheap to construct and apply, J-D can
 *  be very effective.  In the best case, it achieves cubic and quadratic
 *  convergence for symmetric and non-symmetric problems.  This is because
 *  using exact terms, J-D becomes shifted inverse iteration for non-symmetric
 *  problems and Rayleigh quotient iteration for symmetric problems.
 *
 *  Note, SLEPc also has J-D implemented.  We implement it here because it
 *  is relatively simple and it gives an option beyond @ref PowerIteration.
 *  Only the (real part of the) dominant mode is found.
 */

class JacobiDavidson: public EigenSolver
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef EigenSolver                       Base;
  typedef Base::SP_matrix                   SP_matrix;
  typedef Base::SP_solver                   SP_solver;
  typedef Base::SP_vector                   SP_vector;
  typedef Preconditioner::SP_preconditioner SP_pc;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  JacobiDavidson(const double    tol = 1e-6,
                 const int       maxit = 100,
                 const int       subspace_size = 20);

  virtual ~JacobiDavidson(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Sets the operators for the problem.
   *  @param A      left side operator
   *  @param B      optional right side operator (to be inverted)
   *  @param db     optional database for solver and preconditioner options
   */
  void set_operators(SP_matrix A,
                     SP_matrix B = SP_matrix(0),
                     SP_db     db = SP_db(0));

  /**
   *   @brief Sets the preconditioner for the problem
   *   @param P     preconditioner operator (a matrix object)
   */
  void set_preconditioner(SP_matrix P);

  double get_theta() const {return d_theta;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Subspace size
  size_t d_subspace_size;
  /// Preconditioner, i.e. approximate action inv(A-theta*I) * v
  SP_pc d_P;
  /// Current best estimate eigenvalue
  double d_theta;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  void solve_impl(Vector &x, Vector &x0);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

};

/// Default preconditioning process
class JacobiDavidsonDefaultP: public PCShell, public MatrixShell
{

public:

  JacobiDavidsonDefaultP(MatrixBase::SP_matrix        A,
                         MatrixBase::SP_matrix        B,
                         EigenSolver::SP_linearsolver solver,
                         void*                        context)
    : PCShell("jacobi-davidson-default", context)
    , MatrixShell(NULL, A->number_rows(), A->number_columns())
    , d_A(A)
    , d_B(B)
    , d_solver(solver)
    , d_dummy(d_n, 0.0)
  {
    d_solver->set_operators(d_A, d_B);
  }

  /// Perform (A-theta*x)*v or if applicable (inv(B)*A-theta*x)*v
  void multiply(const Vector &x,  Vector &y)
  {
    double theta = ((JacobiDavidson*)PCShell::d_context)->get_theta();
    // dummy <-- A*x
    d_A->multiply(x, d_dummy);
    // y  <-- B\A*x
    if (d_B)
      d_solver->solve(d_dummy, y);
    else
      y.copy(d_dummy);
    // y <-- B\A*x - theta*x = (B\A-I*theta)*x
    y.add_a_times_x(-theta, x);
  }

  void multiply_transpose(const Vector &x,  Vector &y)
  {
    THROW("NOT IMPLEMENTED")
  }

  /// Apply the preconditioner
  void apply(Vector &x,  Vector &y)
  {
    d_solver->solve(x, y);
  }

private:

  /// Right hand side operator
  MatrixBase::SP_matrix d_A;
  /// Left hand side operator (possibly NULL)
  MatrixBase::SP_matrix d_B;
  /// Solver to invert B
  EigenSolver::SP_linearsolver d_solver;
  /// Dummy vector
  Vector d_dummy;
};

} // end namespace callow

//#include "JacobiDavidson.i.hh"

#endif /* callow_JACOBIDAVIDSON_HH_ */
