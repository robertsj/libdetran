//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_SD_1D.i.hh
 * \author robertsj
 * \date   Jun 9, 2012
 * \brief  Equation_SD_1D.i class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SD_1D_I_HH_
#define EQUATION_SD_1D_I_HH_

#include <iostream>

namespace detran
{

inline void Equation_SD_1D::solve(int i,
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
  Require(j == 0);
  Require(k == 0);

  // Compute cell-center angular flux.
  int cell = d_mesh->index(i);
  double coef = 1.0 / (d_material->sigma_t(d_mat_map[cell], d_g) +
                       d_coef_x[i]);
  double psi_center = coef * (source[cell] + d_coef_x[i] * psi_in);

  // Compute outgoing fluxes.
  psi_out = psi_center;

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi) psi[cell] = psi_center;

}

} // end namespace detran

#endif /* EQUATION_SD_1D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SD_1D.i.hh
//---------------------------------------------------------------------------//
