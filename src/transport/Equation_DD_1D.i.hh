//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_DD_1D.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_DD_1D inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_DD_1D_I_HH_
#define detran_EQUATION_DD_1D_I_HH_

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
inline void Equation_DD_1D::setup_angle(const size_t angle)
{
  Require(angle < d_quadrature->number_angles_octant());
  double mu = d_quadrature->mu(0, angle);
  for (size_t i = 0; i < d_mesh->number_cells_x(); ++i)
  {
    d_coef_x[i] = 2.0 * mu / d_mesh->dx(i);
  }
  d_angle = angle;
}

//---------------------------------------------------------------------------//
inline void Equation_DD_1D::solve(const size_t       i,
                                  const size_t       j,
                                  const size_t       k,
                                  moments_type      &source,
                                  face_flux_type    &psi_in,
                                  face_flux_type    &psi_out,
                                  moments_type      &phi,
                                  angular_flux_type &psi)
{
  Require(j == 0);
  Require(k == 0);

  // Compute cell-center angular flux.
  size_t cell = d_mesh->index(i);
  double coef = 1.0 /
                (d_material->sigma_t(d_mat_map[cell], d_g) + d_coef_x[i]);
  double psi_center = coef * (source[cell] + d_coef_x[i] * psi_in);

  // Compute outgoing fluxes.
  psi_out = 2.0*psi_center - psi_in;

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi) psi[cell] = psi_center;

}

} // end namespace detran

#endif /* detran_EQUATION_DD_1D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_1D.i.hh
//---------------------------------------------------------------------------//
