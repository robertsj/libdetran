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
                             SP_db db = SP_db(0))
  {
    // Preconditions
    Require(A);
    d_A = A;
    Ensure(d_A->number_rows() == d_A->number_columns());

    // Default PC settings.
    std::string pc_type = "";
    std::string petsc_pc_type = "ilu";
    int pc_side = LEFT;

    // Extract the PETSc preconditioner.
    PC pc;
    PetscErrorCode ierr;
    ierr = KSPGetPC(d_petsc_solver, &pc);

    // Check database options
    if (db)
    {
      // Get the type
      if(db->check("pc_type"))
        pc_type = db->get<std::string>("pc_type");

      // Check for a callow pc type
      if (pc_type != "petsc_pc")
      {
        if (pc_type == "ilu0")
        {
          d_P = new PCILU0(d_A);
        }
        else if (pc_type == "jacobi")
        {
          d_P = new PCJacobi(d_A);
        }
        // Set callow pc as a shell and set the shell operator
        ierr = PCSetType(pc, PCSHELL);
        d_P->set_petsc_pc(pc);
      }
      // Check for a PETSc pc type
      else
      {
        if(db->check("petsc_pc_type"))
          petsc_pc_type = db->get<std::string>("petsc_pc_type");

        // ILU(k)
        if (petsc_pc_type == "ilu")
        {
          PCSetType(pc, PCILU);
          // Set levels
          int levels = 2;
          if(db->check("petsc_pc_factor_levels"))
            levels = db->get<int>("petsc_pc_factor_levels");
          PCFactorSetLevels(pc, levels);
          // Fill zeros on the diagonal, even if they wouldn't be
          PCFactorSetAllowDiagonalFill(pc);
        }
        // Add more options, like hypre etc.
        else
        {
          THROW("Unsupported PETSc preconditioner: " + petsc_pc_type);
        }
      }
      if(db->check("pc_side"))
        pc_side = db->get<int>("pc_side");
    }
    // If not database, we do not set the preconditioner, letting
    // PETSc do its default.
    {
      PCSetType(pc, PCNONE);
    }

    // Set the operator.  Note, this structure really eliminates
    // the option of using a second matrix for preconditioning.
    // That's fine, though, since any non-A pc's will be shell
    // operations, which are independent of this choice.
    ierr = KSPSetOperators(d_petsc_solver,
                           d_A->petsc_matrix(),
                           d_A->petsc_matrix(),
                           SAME_NONZERO_PATTERN);
    KSPGMRESSetRestart(d_petsc_solver, 20);
    // Set the preconditioner side.
    if (pc_side == Base::LEFT)
      ierr = KSPSetPCSide(d_petsc_solver, PC_LEFT);
    else
      ierr = KSPSetPCSide(d_petsc_solver, PC_RIGHT);

    // Postconditions
    Ensure(!ierr);
  }

  /**
   *  Set the preconditioner.  This allows the client to build, change, etc.
   */
  virtual void set_preconditioner(SP_preconditioner P, const int side = LEFT)
  {
    // Preconditions
    Require(P);
    Require(side == LEFT or side == RIGHT);

    // Keep the pc and its side
    d_P = P;
    d_pc_side = side;

    // Extract the PETSc preconditioner and set P as a shell
    PC pc;
    PetscErrorCode ierr;
    ierr = KSPGetPC(d_petsc_solver, &pc);
    ierr = PCSetType(pc, PCSHELL);
    d_P->set_petsc_pc(pc);

    // Set the preconditioner side.
    if (side == Base::LEFT)
      ierr = KSPSetPCSide(d_petsc_solver, PC_LEFT);
    else
      ierr = KSPSetPCSide(d_petsc_solver, PC_RIGHT);

    // Postconditions
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
