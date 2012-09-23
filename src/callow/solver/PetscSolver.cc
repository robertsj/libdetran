//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PetscSolver.cc
 * \brief  PetscSolver 
 * \author Jeremy Roberts
 * \date   Sep 20, 2012
 */
//---------------------------------------------------------------------------//

#include "callow_config.hh"

#ifdef CALLOW_ENABLE_PETSC

#include "PetscSolver.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR & DESTRUCTOR
//---------------------------------------------------------------------------//

PetscSolver::PetscSolver(const double  atol,
                         const double  rtol,
                         const int     maxit)
  : LinearSolver<PetscScalar>(atol, rtol, maxit, "solver_petsc")
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
}

PetscSolver::~PetscSolver()
{
  KSPDestroy(&d_petsc_solver);
}


} // end namespace callow

#endif

//---------------------------------------------------------------------------//
//              end of file PetscSolver.cc
//---------------------------------------------------------------------------//
