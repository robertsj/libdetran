//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerWG.i.hh
 * \brief  PreconditionerWG.i 
 * \author Jeremy Roberts
 * \date   Aug 5, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONERWG_I_HH_
#define PRECONDITIONERWG_I_HH_

namespace detran
{

// Perform the action y <-- (I + inv(C)*S)*x.  Note the larger comments
// follow the steps given by Larsen in Morel in Nuclear Computational Science.
inline PetscErrorCode PreconditionerWG::apply(Vec x, Vec y)
{
  PetscErrorCode ierr;

  // Get the arrays from x and y, for eventual placement in d_x and d_y.  We
  // do this since d_x an d_y are sized to the diffusion solver, while
  // x and y might have extra unknowns for boundaries.  By simply swapping
  // arrays, we can operate on only the part of x and y needed.  Any remaining
  // part of x (for boundaries) is simply copied to y.
  double *x_a;
  double *y_a;
  ierr = VecGetArray(x, &x_a);
  ierr = VecGetArray(y, &y_a);
  ierr = VecPlaceArray(d_y, y_a);

  // Copy boundary part of x to y.
  int size_full;
  ierr = VecGetSize(x, &size_full);
  int size_moments = d_mesh->number_cells();
  for (int i = size_moments; i < size_full; i++)
    y_a[i] = x_a[i];

  // Temporary vectors
  State::moments_type S_times_x_v(size_moments, 0.0);
  State::moments_type x_v(size_moments, 0.0);

  //-------------------------------------------------------------------------//
  // BEGIN WITH V0.  (Given x, copy into xv for manipulation)
  //-------------------------------------------------------------------------//

  for (int i = 0; i < size_moments; i++)
    x_v[i] = x_a[i];

  //-------------------------------------------------------------------------//
  // OPERATE: V1 <-- S*V0.
  //-------------------------------------------------------------------------//

  d_scattersource->build_within_group_source(d_group, x_v, S_times_x_v);

  // Place S_times_xv array into x.  Note that x has become the
  // right hand side for C*x' = S*x = B.
  ierr = VecPlaceArray(d_x, &S_times_x_v[0]);

  //-------------------------------------------------------------------------//
  // OPERATE: V2 <-- inv(C)*V1 = inv(C)*S*V0.
  //-------------------------------------------------------------------------//

  // Solve(ksp, B, X)
  ierr = KSPSolve(d_solver[d_group], d_x, d_y);
  Insist(!ierr, "Error in KSPSolve.");

  // Reset d_x, and place x's array into d_x
  ierr = VecResetArray(d_x);
  ierr = VecPlaceArray(d_x, x_a);

  //-------------------------------------------------------------------------//
  // OPERATE: V3 <-- V0 + V2 = V0 + inv(C)*S*V0 = (I + inv(C)*S)*V0
  //-------------------------------------------------------------------------//

  // Do y <-- x + y = x + inv(C)*S*x
  ierr = VecAXPY(d_y, 1.0, d_x);

  // Reset the Vec's
  VecResetArray(d_y);
  VecResetArray(d_x);
  VecRestoreArray(y, &y_a);
  VecRestoreArray(x, &x_a);

  return ierr;
}

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

inline PetscErrorCode apply_inv_P(PC wg_pc, Vec x, Vec y)
{
//  // Get the context and cast as InnerGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = PCShellGetContext(wg_pc, &ctx); CHKERRQ(ierr);
  detran::PreconditionerWG *tmp =
    (detran::PreconditionerWG*) ctx;
  // Call the actual apply operator.
  return tmp->apply(x, y);
}

} // end namespace detran

#endif // PRECONDITIONERWG_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file PreconditionerWG.i.hh
//---------------------------------------------------------------------------//
