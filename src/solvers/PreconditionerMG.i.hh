//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerMG.i.hh
 * \brief  PreconditionerMG.i 
 * \author Jeremy Roberts
 * \date   Aug 6, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PRECONDITIONERWG_I_HH_
#define PRECONDITIONERWG_I_HH_

namespace detran
{

// Perform the action y <-- (I + inv(C)*S)*x.  Note the larger comments
// follow the steps given by Larsen in Morel in Nuclear Computational Science.
inline PetscErrorCode PreconditionerMG::apply(Vec x, Vec y)
{
  PetscErrorCode ierr;

  // Get the arrays from x and y, for eventual placement in d_x and d_y.  We
  // do this since d_x an d_y are sized to the diffusion solver, while
  // x and y might have extra unknowns for boundaries.  By simply swapping
  // arrays, we can operate on only the part of x and y needed.  Any remaining
  // part of x (for boundaries) is simply copied to y.
  Vec z;
  VecDuplicate(x, &z);
  double *x_a;
  double *y_a;
  double *z_a;
  ierr = VecGetArray(x, &x_a);
  ierr = VecGetArray(y, &y_a);
  ierr = VecGetArray(z, &z_a);

  // build_total_group_source(int g, int g_cutoff,
  // const State::vec_moments_type &phi,
  // moments_type &source)

  int ng = d_material->number_groups();

  // Copy x_a to the vec_moments_type.  Only the upscatter block is used.
  State::vec_moments_type
    phi(ng, State::moments_type(d_moments_size_group, 0.0));
  for (int g = d_upscatter_cutoff; g < ng; g++)
  {
    for (int i = 0; i < d_moments_size_group; i++)
    {
      phi[g][i] = x_a[(g - d_upscatter_cutoff) * d_moments_size_group + i];
    }
  }

  // Create the total group source, and copy into x_a.
  State::moments_type source(d_moments_size_group, 0.0);
  for (int g = d_upscatter_cutoff; g < ng; g++)
  {
    d_scattersource->build_total_group_source(g, d_upscatter_cutoff, phi, source);
    for (int i = 0; i < d_moments_size_group; i++)
    {
      z_a[(g - d_upscatter_cutoff) * d_moments_size_group + i] = source[i];
    }
  }
  // x_a now has S_times_x

  //-------------------------------------------------------------------------//
  // OPERATE: V2 <-- inv(C)*V1 = inv(C)*S*V0.
  //-------------------------------------------------------------------------//

  // Solve(ksp, B, X)
  ierr = KSPSolve(d_solver, z, y);
  Insist(!ierr, "Error in KSPSolve.");

  //-------------------------------------------------------------------------//
  // OPERATE: V3 <-- V0 + V2 = V0 + inv(C)*S*V0 = (I + inv(C)*S)*V0
  //-------------------------------------------------------------------------//

  // Do y <-- x + y = x + inv(C)*S*x
  ierr = VecAXPY(y, 1.0, x);

  // Reset the Vec's
  VecRestoreArray(z, &z_a);
  VecRestoreArray(y, &y_a);
  VecRestoreArray(x, &x_a);

  VecDestroy(&z);
  return ierr;
}

//---------------------------------------------------------------------------//
// EXTERNAL WRAPPER FUNCTIONS
//---------------------------------------------------------------------------//

inline PetscErrorCode apply_inv_P_MG(PC wg_pc, Vec x, Vec y)
{
//  // Get the context and cast as InnerGMRES pointer.
  PetscErrorCode ierr;
  void *ctx;
  ierr = PCShellGetContext(wg_pc, &ctx); CHKERRQ(ierr);
  detran::PreconditionerMG *tmp =
    (detran::PreconditionerMG*) ctx;
  // Call the actual apply operator.
  return tmp->apply(x, y);
}

} // end namespace detran

#endif // PRECONDITIONERWG_I_HH_

//---------------------------------------------------------------------------//
//              end of file PreconditionerMG.i.hh
//---------------------------------------------------------------------------//
