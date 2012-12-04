//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_1D.i.hh
 *  @author robertsj
 *  @date   Jun 9, 2012
 *  @brief  Equation_SC_1D.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SC_1D_I_HH_
#define EQUATION_SC_1D_I_HH_

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
inline void Equation_SC_1D::setup_angle(const size_t angle)
{
  Require(angle < d_quadrature->number_angles_octant());
  d_angle = angle;
  d_mu = d_quadrature->mu(0, d_angle);
}

//---------------------------------------------------------------------------//
inline void Equation_SC_1D::solve(const size_t i,
                                  const size_t j,
                                  const size_t k,
                                  moments_type &source,
                                  face_flux_type &psi_in,
                                  face_flux_type &psi_out,
                                  moments_type &phi,
                                  angular_flux_type &psi)
{
  // Preconditions.  (The client *must* set group and angles.)
  Require(j == 0);
  Require(k == 0);
  Require(d_mu > 0.0);

  // Compute cell-center angular flux.
  double sigma = d_material->sigma_t(d_mat_map[i], d_g);
  double tau   = sigma * d_mesh->dx(i) / d_mu;
  double A     = std::exp(-tau);
  double q     = source[i];

  // Cell average flux
  double psi_avg = psi_in * (1.0 - A) / tau +
                   q * (sigma * d_mesh->dx(i) + d_mu * (A - 1.0)) /
                       (sigma * sigma * d_mesh->dx(i));

  // Compute outgoing fluxes.
  psi_out = A * psi_in + q * (1.0 - A) / sigma;

  // Compute flux moments.
  phi[i] += d_quadrature->weight(d_angle) * psi_avg;

  // Store angular flux if needed.
  if (d_update_psi) psi[i] = psi_avg;

}

} // end namespace detran

#endif /* EQUATION_SC_1D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_1D.i.hh
//---------------------------------------------------------------------------//
