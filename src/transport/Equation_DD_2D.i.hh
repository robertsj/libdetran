//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_2D.i.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_2D inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_DD_2D_I_HH_
#define EQUATION_DD_2D_I_HH_

namespace detran
{

template <>
void Equation_DD_2D::solve(int g,
                   int i,
                   int j,
                   int k,
                   face_flux_type psi_in,
                   double source,
                   face_flux_type psi_out,
                   moments_type &phi,
                   angular_flux_type &psi)
{
  // Compute cell-center angular flux.
  int cell = d_mesh->index(i, j);
  double coef = 1.0 / (d_material->sigma_t(d_mat_map[cell], g) +
                       d_coef_x[i] + d_coef_y[j]);
  double psi_center = coef * (source + d_coef_x[i] * psi_in[1] +
                                       d_coef_y[j] * psi_in[0] );

  // Compute outgoing fluxes.
  psi_out[0] = 2.0*psi_center - psi_in[0];
  psi_out[1] = 2.0*psi_center - psi_in[1];

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle);

  // Store angular flux if needed.
  if (d_store_psi)
  {
    psi[cell] = psi_center;
  }
}

} // end namespace detran

#endif /* EQUATION_DD_2D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_2D.i.hh
//---------------------------------------------------------------------------//
