//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PetscSolver.cc
 *  @brief PetscSolver member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "callow_config.hh"

#ifdef CALLOW_ENABLE_PETSC

#include "PetscSolver.hh"

namespace callow
{

//----------------------------------------------------------------------------//
PetscSolver::PetscSolver(const double  atol,
                         const double  rtol,
                         const int     maxit)
  : LinearSolver(atol, rtol, maxit, "petsc")
{

  // create the petsc object
  PetscErrorCode ierr = KSPCreate(PETSC_COMM_SELF, &d_petsc_solver);
  Insist(!ierr, "Error creating KSP object.");

  // set the convergence criteria (third param is divergence tolerace; skip)
  ierr = KSPSetTolerances(d_petsc_solver, rtol, atol, PETSC_DEFAULT, maxit);
  Insist(!ierr, "Error setting KSP tolerances.");

  // Allow for command line flags.
  ierr = KSPSetFromOptions(d_petsc_solver);
  Insist(!ierr, "Error setting KSP from options.");

  /// set the monitor so we can extract residuals on-the-fly
  ierr = KSPMonitorSet(d_petsc_solver, petsc_ksp_monitor, this, PETSC_NULL);
  Insist(!ierr, "Error setting KSP monitor.");

  ierr = KSPSetInitialGuessNonzero(d_petsc_solver, PETSC_TRUE);
  Insist(!ierr, "Error setting KSP nonzeros initial guess.");
}

//----------------------------------------------------------------------------//
PetscSolver::~PetscSolver()
{
  KSPDestroy(&d_petsc_solver);
}

//----------------------------------------------------------------------------//
void PetscSolver::set_operators(SP_matrix A, SP_db db)
{
  // Preconditions
  Require(A);
  d_A = A;
  Ensure(d_A->number_rows() == d_A->number_columns());

  // Set the db, if present.  Otherwise, d_db is unchanged, which
  // lets us set new operators but maintain old parameters.
  if (db) d_db = db;

  // Default PC settings.
  std::string pc_type = "";
  std::string petsc_pc_type = "ilu";
  int pc_side = RIGHT;

  // Extract the PETSc preconditioner.
  PC pc;
  PetscErrorCode ierr;
  ierr = KSPGetPC(d_petsc_solver, &pc);

  // Check database options
  if (d_db)
  {
    // Get the type
    if (d_db->check("pc_type")) pc_type = d_db->get<std::string>("pc_type");

    // Check for a callow pc type
    if (pc_type != "petsc_pc")
    {
      if (pc_type == "ilu0")
        d_P = new PCILU0(d_A);
      else if (pc_type == "jacobi")
        d_P = new PCJacobi(d_A);
      // Set callow pc as a shell and set the shell operator
      if (d_P)
      {
        ierr = PCSetType(pc, PCSHELL);
        d_P->set_petsc_pc(pc);
      }
    }
    // Check for a PETSc pc type
    else
    {
      if (d_db->check("petsc_pc_type"))
        petsc_pc_type = d_db->get<std::string>("petsc_pc_type");

      // ILU(k)
      if (petsc_pc_type == "ilu")
      {
        PCSetType(pc, PCILU);
        // Set levels
        int levels = 2;
        if (d_db->check("petsc_pc_factor_levels"))
          levels = d_db->get<int>("petsc_pc_factor_levels");
        PCFactorSetLevels(pc, levels);
        // Fill zeros on the diagonal, even if they wouldn't be
        PCFactorSetAllowDiagonalFill(pc);
      }
      // LU
      else if (petsc_pc_type == "lu")
      {
        PCSetType(pc, PCLU);
        PCFactorSetReuseOrdering(pc, PETSC_TRUE);
        PCFactorSetReuseFill(pc, PETSC_TRUE);
//        PCFactorSetMatOrderingType(pc, MATORDERINGND);
//        PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU);
      }
      // Add more options, like hypre etc.
      else
      {
        //PCSetType(pc, petsc_pc_type.c_str());
        PCSetFromOptions(pc);
        //THROW("Unsupported PETSc preconditioner: " + petsc_pc_type);
      }
    }
    if (d_db->check("pc_side")) pc_side = d_db->get<int>("pc_side");
  }
  // If not database, we do not set the preconditioner, letting
  // PETSc do its default.
  else
  {
    //std::cout << " WHY NO PETSC HAS PC? " << std::endl;
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

  KSPGMRESSetRestart(d_petsc_solver, 30);
  // Set the preconditioner side.
  if (pc_side == Base::LEFT) ierr = KSPSetPCSide(d_petsc_solver, PC_LEFT);
  else
    ierr = KSPSetPCSide(d_petsc_solver, PC_RIGHT);

  // Postconditions
  Ensure(!ierr);
}

//----------------------------------------------------------------------------//
void PetscSolver::set_preconditioner(SP_preconditioner P,
                                     const int         side)
{
  Require(P);
  Require(side == LEFT || side == RIGHT);

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

  Ensure(!ierr);
}

} // end namespace callow

#endif

//----------------------------------------------------------------------------//
//              end of file PetscSolver.cc
//----------------------------------------------------------------------------//
