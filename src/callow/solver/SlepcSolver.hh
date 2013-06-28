//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SlepcSolver.hh
 *  @brief SlepcSolver class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_SLEPCSOLVER_HH_
#define callow_SLEPCSOLVER_HH_

#include "EigenSolver.hh"

namespace callow
{

/**
 *  @class SlepcSolver
 *  @brief Solve the eigenvalue problem with SLEPc
 */
class SlepcSolver: public EigenSolver
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef EigenSolver                       Base;
  typedef Base::SP_matrix                   SP_matrix;
  typedef Base::SP_solver                   SP_solver;
  typedef Base::SP_vector                   SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  SlepcSolver(const std::string &epstype,
              const double       tol,
              const int          maxit,
              const int          number_values);

  virtual ~SlepcSolver();

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
  void set_operators(SP_matrix A,
                     SP_matrix B  = SP_matrix(0),
                     SP_db     db = SP_db(0));


  /// Set the preconditioner for a generalized eigenvalue problem
  void set_preconditioner(SP_preconditioner pc,
                          const int         side = LinearSolver::LEFT);

  /// Set the preconditioner matrix for a generalized eigenvalue problem
  void set_preconditioner_matrix(SP_matrix P,
                                 const int side = LinearSolver::LEFT);

  /// Set the EPS type
  void set_eps_type(const std::string &eps_type);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  // expose base members
  using Base::d_tolerance;
  using Base::d_maximum_iterations;
  using Base::d_name;
  using Base::d_residual_norm;
  using Base::d_number_iterations;
  using Base::d_A;
  using Base::d_B;
  using Base::d_solver;
  using Base::d_monitor_level;
  using Base::d_lambda;

  /// SLEPc solver type
  std::string d_eps_type;
  /// SLEPc solver object
  EPS d_slepc_solver;
  /// Optional preconditioner
  SP_preconditioner d_P;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EIGENSOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  virtual void solve_impl(Vector &x, Vector &x0);

};

} // end namespace callow

#endif /* callow_SLEPCSOLVER_HH_ */

//----------------------------------------------------------------------------//
//              end of file SlepcSolver.hh
//----------------------------------------------------------------------------//
