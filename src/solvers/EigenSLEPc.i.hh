//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   EigenSLEPc.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  EigenSLEPc inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EIGENSLEPC_I_HH_
#define EIGENSLEPC_I_HH_

// Utilities
#include "Warning.hh"

// System
#include <algorithm>
#include <cstdio>
#include <iostream>

namespace detran
{

template <class D>
void EigenSLEPc<D>::solve()
{

  using std::cout;
  using std::endl;

  cout << "Starting EIGEN." << endl;

  // Initialize the fission density
  b_fissionsource->initialize();

  // Temporaries
  double lambda;
  double lambda_imag;
  Vec    J_real;
  Vec    J_imag;
  int    ierr;

  ierr = EPSSetTolerances(d_solver, b_tolerance, b_max_iters);
  Insist(!ierr, "Error setting EPS tolerances.");

  ierr = EPSSolve(d_solver);
  Insist(!ierr, "Error solving eigenvalue problem.");

  int numit = 0;
  ierr = EPSGetIterationNumber(d_solver, &numit);
  Insist(!ierr, "Error getting iteration count.");
  std::cout << "EIGEN: Number of iterations =" << numit << std::endl;

  // Get the dominant mode.
  MatGetVecs(d_operator, PETSC_NULL, &J_imag);
  MatGetVecs(d_operator, PETSC_NULL, &J_real);
  ierr = EPSGetEigenpair(d_solver, 0, &lambda, &lambda_imag, J_real, J_imag);
  Insist(!ierr, "Error getting eigenpair.");

  // Scale the result by its sum.  This points it in the positive
  // direction and gives the L1 normalization we're using in PI.
  double sum;
  VecSum(J_real, &sum);
  VecScale(J_real, 1.0/sum);

  // Copy the result into the fission source.
  VecPlaceArray(J_imag, &(b_fissionsource->density()[0]));
  VecCopy(J_real, J_imag);
  VecResetArray(J_imag);

  // To retrieve the correct flux moments, we need to do
  // one more solve with this new density.  This could
  // be a switched feature.
  b_mg_solver->solve();

  // Free temporary
  VecDestroy(&J_imag);
  VecDestroy(&J_real);

  b_state->set_eigenvalue(lambda);
  std::cout << "EIGEN done. keff = " << lambda <<  std::endl;
  std::cout << "Number MG solves = " << d_mg_solves <<  std::endl;
}

template <class D>
inline PetscErrorCode EigenSLEPc<D>::apply_eigen(Mat A, Vec X, Vec Y)
{
  // Increment the counter.
  d_mg_solves++;

  // Temporarily swap Y's internal with the density.
  VecPlaceArray(Y, &(b_fissionsource->density()[0]));

  // Copy X into Y.  Now density has the right values.
  VecCopy(X, Y);

  // Reset Y.
  VecResetArray(Y);

  // Solve the multigroup equations.
  b_mg_solver->solve();

  // Update the density.
  b_fissionsource->update();

  // Temporarily swap Y's internal with the density.
  VecPlaceArray(X, &(b_fissionsource->density()[0]));

  // Copy X into Y.  Now Y has the right values.
  VecCopy(X, Y);

  // Reset X.
  VecResetArray(X);

  return 0;
}

// Explicit instantiations
template class EigenSLEPc<_1D>;
template class EigenSLEPc<_2D>;
template class EigenSLEPc<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

inline PetscErrorCode apply_eigen_1D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as EigenSLEPc pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::EigenSLEPc<detran::_1D> *inner =
    (detran::EigenSLEPc<detran::_1D>*) ctx;
  // Call the actual apply operator.
  return inner->apply_eigen(A, x, y);
}

inline PetscErrorCode apply_eigen_2D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as EigenSLEPc pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::EigenSLEPc<detran::_2D> *inner =
    (detran::EigenSLEPc<detran::_2D>*) ctx;
  // Call the actual apply operator.
  return inner->apply_eigen(A, x, y);
}

inline PetscErrorCode apply_eigen_3D(Mat A, Vec x, Vec y)
{
  // Get the context and cast as EigenSLEPc pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = MatShellGetContext(A, &ctx); CHKERRQ(ierr);
  detran::EigenSLEPc<detran::_3D> *inner =
    (detran::EigenSLEPc<detran::_3D>*) ctx;
  // Call the actual apply operator.
  return inner->apply_eigen(A, x, y);
}

#endif /* EIGENSLEPC_I_HH_ */
