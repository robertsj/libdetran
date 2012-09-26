//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SlepcSolver.i.hh
 *  @author robertsj
 *  @date   Sep 25, 2012
 *  @brief  SlepcSolver inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef callow_SLEPCSOLVER_I_HH_
#define callow_SLEPCSOLVER_I_HH_

namespace callow
{

inline void SlepcSolver::
solve_impl(Vector<PetscScalar> &x, Vector<PetscScalar> &x0)
{

  // Temporaries
  double lambda;
  double lambda_imag;

  // Solve the eigenproblem
  PetscErrorCode ierr = EPSSolve(d_slepc_solver);
  Insist(!ierr, "Error solving eigenvalue problem.");

  // Extract the number of iterations
  ierr = EPSGetIterationNumber(d_slepc_solver, &d_number_iterations);
  Insist(!ierr, "Error getting iteration count.");

  // Get the dominant mode.
  ierr = EPSGetEigenpair(d_slepc_solver, 0, &lambda, &lambda_imag,
                         x.petsc_vector(), x0.petsc_vector());
  Insist(!ierr, "Error getting eigenpair.");

  // Scale the result by its sum.  This points it in the positive
  // direction and gives the L1 normalization we're using in PI.
  PetscScalar sign = 1;
  if (x[0] < 0) sign = -1;
  x.scale(sign / x.norm(L1));

  // Store the eigenvalue
  d_lambda = lambda;

}


} // end namespace callow

#endif /* callow_SLEPCSOLVER_I_HH_ */
