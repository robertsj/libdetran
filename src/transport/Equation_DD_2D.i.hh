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

#include <iostream>

namespace detran
{

inline void Equation_DD_2D::solve(int i,
                                  int j,
                                  int k,
                                  moments_type &source,
                                  face_flux_type &psi_in,
                                  face_flux_type &psi_out,
                                  moments_type &phi,
                                  angular_flux_type &psi)
{
  // Preconditions.  (The client *must* set group and angles.)
  Require(d_g >= 0);
  Require(d_angle >= 0);
  Require(d_octant >= 0);
  Require(k == 0);

  // Compute cell-center angular flux.
  int cell = d_mesh->index(i, j);
  double coef = 1.0 / (d_material->sigma_t(d_mat_map[cell], d_g) +
                       d_coef_x[i] + d_coef_y[j]);
  double psi_center = coef * (source[cell] + d_coef_x[i] * psi_in[Mesh::VERT] +
                                             d_coef_y[j] * psi_in[Mesh::HORZ] );

  // Compute outgoing fluxes.
  psi_out[Mesh::HORZ] = 2.0*psi_center - psi_in[Mesh::HORZ];
  psi_out[Mesh::VERT] = 2.0*psi_center - psi_in[Mesh::VERT];

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi)
  {
    psi[cell] = psi_center;
  }
}

} // end namespace detran

#endif /* EQUATION_DD_2D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_2D.i.hh
//---------------------------------------------------------------------------//
