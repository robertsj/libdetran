//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_3D.i.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_3D inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_DD_3D_I_HH_
#define EQUATION_DD_3D_I_HH_

#include <iostream>

namespace detran
{

inline void Equation_DD_3D::solve(int i,
                                  int j,
                                  int k,
                                  moments_type &source,
                                  face_flux_type &psi_in,
                                  face_flux_type &psi_out,
                                  moments_type &phi,
                                  angular_flux_type &psi)
{
  using std::cout;
  using std::endl;

  // Preconditions.  (The client *must* set group and angles.)
  Require(d_g >= 0);
  Require(d_angle >= 0);
  Require(d_octant >= 0);

  // Compute cell-center angular flux.
  int cell = d_mesh->index(i, j, k);
  double coef = 1.0 / (d_material->sigma_t(d_mat_map[cell], d_g) +
                       d_coef_x[i] + d_coef_y[j] + d_coef_z[k]);
  double psi_center = coef * (source[cell] + d_coef_x[i] * psi_in[Mesh::YZ] +
                                             d_coef_y[j] * psi_in[Mesh::XZ] +
                                             d_coef_z[k] * psi_in[Mesh::XY]);

  // Compute outgoing fluxes.
  double two_psi_center = 2.0 * psi_center;
  psi_out[0] = two_psi_center - psi_in[0];
  psi_out[1] = two_psi_center - psi_in[1];
  psi_out[2] = two_psi_center - psi_in[2];

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

//  cout << "------------------------------------" << endl;
//  cout << " [i, j, k] = " << " [" << i << "," << j << "," << k << "]" << endl;
//  cout << " cell = " << cell << endl;
//  cout << " source = " << source[cell] << endl;
//  cout << "------------------------------------" << endl;

  // Store angular flux if needed.
  if (d_update_psi)
  {
    psi[cell] = psi_center;
  }
}

} // end namespace detran

#endif /* EQUATION_DD_3D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_3D.i.hh
//---------------------------------------------------------------------------//
