//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PetscSolver.hh
 *  @author robertsj
 *  @date   Sep 20, 2012
 *  @brief  PetscSolver class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_PETSCSOLVER_HH_
#define callow_PETSCSOLVER_HH_

#ifdef CALLOW_ENABLE_PETSC

#include "LinearSolver.hh"

namespace callow
{

/**
 *  @class PetscSolver
 *  @brief Uses PETSc to solve a system
 */
class PetscSolver: public LinearSolver
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef LinearSolver                          Base;
  typedef Base::SP_solver                       SP_solver;
  typedef MatrixBase ::SP_matrix                SP_matrix;
  typedef Preconditioner::SP_preconditioner     SP_preconditioner;
  typedef Vector::SP_vector                     SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  PetscSolver(const double atol, const double rtol, const int maxit);

  virtual ~PetscSolver();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  This overloads the default implementation so that we can extract
   *  the PETSc PC object and set it in our PC object if present.
   *
   */
  virtual void set_operators(SP_matrix A,
                             SP_preconditioner P = SP_preconditioner(0),
                             const int side = Base::LEFT)
  {
    Require(A);
    d_A = A;
    Ensure(d_A->number_rows() == d_A->number_columns());
    PetscErrorCode ierr;

    // extract the preconditioner
    PC pc;
    ierr = KSPGetPC(d_petsc_solver, &pc);

    // if a preconditioner object is present, set it
    if (P)
    {
      d_P  = P;
      ierr = PCSetType(pc, PCSHELL);
      // internally, this tells petsc the shell operator
      d_P->set_petsc_pc(pc);
    }
    // else do something else
    {
      ierr = PCSetType(pc, PCILU);
      ierr = PCFactorSetLevels(pc, 2);
    }

    // Set the operator.
    ierr = KSPSetOperators(d_petsc_solver,
                           d_A->petsc_matrix(),
                           d_A->petsc_matrix(),
                           SAME_NONZERO_PATTERN);

    if (side == Base::LEFT)
      ierr = KSPSetPCSide(d_petsc_solver, PC_LEFT);
    else if (side == Base::RIGHT)
      ierr = KSPSetPCSide(d_petsc_solver, PC_RIGHT);
    Ensure(!ierr);
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

  // petsc solver type
  KSP d_petsc_solver;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL LINEAR SOLVERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /**
   *  \param b  right hand side
   *  \param x  unknown vector
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

#endif // CALLOW_ENANLE_PETSC

#endif /* callow_PETSCSOLVER_HH_ */
