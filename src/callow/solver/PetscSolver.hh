//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   PetscSolver.hh
 *  @brief  PetscSolver class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_PETSCSOLVER_HH_
#define callow_PETSCSOLVER_HH_

#ifdef DETRAN_ENABLE_PETSC

#include "LinearSolver.hh"
// preconditioners
#include "callow/preconditioner/PCILU0.hh"
#include "callow/preconditioner/PCJacobi.hh"

namespace callow
{

/**
 *  @class PetscSolver
 *  @brief Uses PETSc to solve a system
 */
class PetscSolver: public LinearSolver
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef LinearSolver                          Base;
  typedef Base::SP_solver                       SP_solver;
  typedef MatrixBase ::SP_matrix                SP_matrix;
  typedef Preconditioner::SP_preconditioner     SP_preconditioner;
  typedef Vector::SP_vector                     SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  PetscSolver(const double atol, const double rtol, const int maxit);

  virtual ~PetscSolver();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  This overloads the default implementation so that we can extract
   *  the PETSc PC object and set it in our PC object if present.
   *
   */
  void set_operators(SP_matrix A, SP_db db = SP_db(0));

  /**
   *  Set the preconditioner.  This allows the client to build, change, etc.
   */
  void set_preconditioner(SP_preconditioner P, const int side = LEFT);

  /// Get the KSP object
  KSP petsc_solver(){return d_petsc_solver;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  // expose base class members
  using LinearSolver::d_absolute_tolerance;
  using LinearSolver::d_relative_tolerance;
  using LinearSolver::d_maximum_iterations;
  using LinearSolver::d_residual;
  using LinearSolver::d_number_iterations;
  using LinearSolver::d_A;
  using LinearSolver::d_P;

  // petsc solver type
  KSP d_petsc_solver;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /**
   *  @param b  right hand side
   *  @param x  unknown vector
   */
  void solve_impl(const Vector &b, Vector &x);

  /// let the monitor wrapper call our monitor
  friend PetscErrorCode
  petsc_ksp_monitor(KSP ksp, PetscInt it, PetscReal rnorm, void* ctx);

};

/// Monitor the solution
PetscErrorCode petsc_ksp_monitor(KSP ksp, PetscInt it,
                                 PetscReal rnorm, void* ctx);

} // end namespace callow

// Inline member definitions
#include "PetscSolver.i.hh"

#endif // DETRAN_ENABLE_PETSC

#endif /* callow_PETSCSOLVER_HH_ */

//----------------------------------------------------------------------------//
//              end of file PetscSolver.hh
//----------------------------------------------------------------------------//
