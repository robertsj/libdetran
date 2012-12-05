//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PetscSolver.i.hh
 *  @author robertsj
 *  @date   Sep 20, 2012
 *  @brief  PetscSolver inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_PETSCSOLVER_I_HH_
#define callow_PETSCSOLVER_I_HH_

#include "callow/matrix/Matrix.hh"

namespace callow
{

//---------------------------------------------------------------------------//
// SOLVE
//---------------------------------------------------------------------------//

inline void PetscSolver::
solve_impl(const Vector &b, Vector &x)
{
  PetscErrorCode ierr;

  // solve the problem
  ierr = KSPSolve(d_petsc_solver,
                  const_cast<Vector*>(&b)->petsc_vector(),
                  x.petsc_vector());

  Insist(!ierr, "Error in KSPSolve.");
}

inline PetscErrorCode
petsc_ksp_monitor(KSP ksp, PetscInt it, PetscReal rnorm, void* ctx)
{
  PetscSolver* solver = (PetscSolver*)ctx;
  solver->monitor(it, rnorm); // note, petsc is in charge of terminating
  return 0;
}

} // end namespace callow

#endif /* callow_PETSCSOLVER_I_HH_ */
