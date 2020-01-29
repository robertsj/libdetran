//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Davidson.hh
 *  @brief Davidson class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_DAVIDSON_HH_
#define callow_DAVIDSON_HH_

#include "EigenSolver.hh"
#include "callow/preconditioner/PCShell.hh"
#include "callow/matrix/MatrixShell.hh"

namespace callow
{

/**
 *  @class Davidson
 *  @brief Solves generalized eigenproblems via generalized Davidson's method
 *
 *  Consider the generalized eigenvalue problem
 *  @f[
 *       \mathbf{A}x = \lambda  \mathbf{B}x \, .
 *  @f]
 *  The basic idea of the Davidson's method is
 *  that at each iteration, the problem is projected into
 *  a subspace, and the resulting Ritz pair is used to estimate the
 *  eigenvalue and eigenvector of the original system.  A correction
 *  to the eigenvector is produced that is (approximately)
 *  orthogonal to the residual.  In other words, for the updated
 *  eigenvector @f$ u @f$ and eigenvalue @f$ \lambda @f$, we solve
 *  @f[
 *     (\mathbf{A} - \lambda \mathbf{B}) v = \mathbf{R} v \approx r \, ,
 *  @f]
 *  or alternatively,
 *  @f[
 *     \mathbf{P} v = r \, ,
 *  @f]
 *  where @f$P \approx R @f$ is called the
 *  <em> preconditioner </em>.
 *  The result @f$ v @f$ is made orthogonal to the current subspace and
 *  is then used to extend subspace.  While not discussed in depth
 *  here, projection methods in general can be much more efficient
 *  if restarted.  We simply take the last resulting eigenpair.
 *
 *  SLEPc also offers an implementation of generalized Davidson
 *  that may be more efficient but cannot handle arbitrary user-defined
 *  preconditioners in the way needed for Detran.
 */
class Davidson: public EigenSolver
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

  Davidson(const double  tol = 1e-6,
           const int     maxit = 100,
           const int     subspace_size = 20);

  virtual ~Davidson(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Sets the operators for the problem.
   *  @param A      left side operator
   *  @param B      right side operator
   *  @param db     optional database for solver and preconditioner options
   */
  void set_operators(SP_matrix A,
                     SP_matrix B  = SP_matrix(0),
                     SP_db     db = SP_db(0));

  /**
   *   @brief Sets the preconditioner for the problem
   *   @param P     preconditioner operator
   *   @param size  preconditioner side (unused for Davidson)
   */
  void set_preconditioner(SP_pc     P,
                          const int side = LinearSolver::LEFT);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Subspace size
  size_t d_subspace_size;
  /// Preconditioner, i.e. approximate action inv(A-theta*B) * v
  SP_pc d_P;
  /// Residual operator
  MatrixBase::SP_matrix d_A_minus_ritz_times_B;
  /// Projected system solver
  SP_solver d_projected_solver;

  //-------------------------------------------------------------0------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  void solve_impl(Vector &x, Vector &x0);

};

/// Wraps the action of r = (A-lambda*B)*x
class DavidsonResidual: public MatrixShell
{

public:

  DavidsonResidual(MatrixBase::SP_matrix A,
                   MatrixBase::SP_matrix B,
                   void*                 context)
    : MatrixShell(context, A->number_rows(), A->number_columns())
    , d_A(A)
    , d_B(B)
    , d_dummy(d_n, 0.0)
  {
    Require(d_A);
    Require(d_B);
  }

  void multiply(const Vector &x,  Vector &y)
  {
    double ritz = ((Davidson*)MatrixShell::d_context)->eigenvalue();
    // dummy <-- A*x
    d_A->multiply(x, y);
    // y  <-- B*x
    d_B->multiply(x, d_dummy);
    // y <-- A*x - ritz*B*x
    y.add_a_times_x(-ritz, d_dummy);
  }

  void multiply_transpose(const Vector &x,  Vector &y)
  {
    THROW("NOT IMPLEMENTED")
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Right hand side operator
  MatrixBase::SP_matrix d_A;
  /// Left hand side operator
  MatrixBase::SP_matrix d_B;
  /// Dummy vector
  Vector d_dummy;

};

/**
 *  @class DavidsonDefaultP
 *  @brief Provides default preconditioning process for nonlinear Arnoldi
 *
 *  The correction is defined by the approximate solve
 *  @f[
 *     \mathbf{A}(\lambda) v \approx r \, .
 *  @f]
 *  Here, a very coarse GMRES solve is performed, using a large tolerance
 *  and low iteration limit.
 */
class DavidsonDefaultP: public PCShell
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Davidson::SP_db        SP_db;
  typedef MatrixBase::SP_matrix  SP_matrix;

  /**
   *  @brief Constructor
   *  @param    A         left hand operator
   *  @param    B         right hand operator
   *  @param    context   user context
   */
  DavidsonDefaultP(MatrixBase::SP_matrix  A_minus_ritz_times_B,
                   void*                  context,
                   SP_db                  db = SP_db(0))
    : PCShell("davidson-default", context)
  {
    if (!db)
    {
      db = detran_utilities::InputDB::Create();
      db->put<std::string>("linear_solver_type",  "gmres");
      db->put<double>("linear_solver_atol",       0.0);
      db->put<double>("linear_solver_rtol",       0.1);
      db->put<int>("linear_solver_maxit",         10);
      db->put<int>("linear_solver_monitor_level", 0);

    }
    d_solver = LinearSolverCreator::Create(db);
    d_solver->set_operators(A_minus_ritz_times_B);
  }

  /// Apply the preconditioner
  void apply(Vector &x,  Vector &y)
  {
    //x.display("X solve");
    // give a 0 initial guess so that (A-e*B)x != y
    y.scale(0.0);
    d_solver->solve(x, y);
    //y.display("Y solve");
  }

private:

  /// Solver for approximating inv(A-e*B)
  EigenSolver::SP_linearsolver d_solver;

};

} // end namespace callow

#endif /* callow_DAVIDSON_HH_ */

//----------------------------------------------------------------------------//
//              end of file Davidson.hh
//----------------------------------------------------------------------------//
