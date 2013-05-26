//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_2D.i.hh
 *  @author robertsj
 *  @date   Apr 10, 2012
 *  @brief  Equation_SC_2D inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_SC_2D_I_HH_
#define detran_EQUATION_SC_2D_I_HH_

// System
#include <cmath>
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
inline void Equation_SC_2D::solve(const size_t       i,
                                  const size_t       j,
                                  const size_t       k,
                                  moments_type      &source,
                                  face_flux_type    &psi_in,
                                  face_flux_type    &psi_out,
                                  moments_type      &phi,
                                  angular_flux_type &psi)
{
  // Preconditions.  (The client *must* set group and angles.)
  Require(i < d_mesh->number_cells_x());
  Require(j < d_mesh->number_cells_y());
  Require(k == 0);

  typedef detran_geometry::Mesh Mesh;

  int cell = d_mesh->index(i, j);
  double sigma = d_material->sigma_t(d_mat_map[cell], d_g);
  double Q = source[cell] / sigma;
  double alpha = sigma * d_alpha[i];
  double beta = sigma * d_beta[j];
  double rho = alpha / beta; // \todo sigma cancels? precompute alpha/beta?
  double psi_in_V_minus_Q = psi_in[Mesh::VERT] - Q;
  double psi_in_H_minus_Q = psi_in[Mesh::HORZ] - Q;

  // Compute outgoing face fluxes.
  if (rho <= 1.0)
  {
    double expf = exp_appx(-alpha);
    double one_m_exp_alpha = (1.0 - expf) / alpha;
    psi_out[Mesh::VERT] = Q + psi_in_V_minus_Q * (1.0 - rho) * expf
                            + psi_in_H_minus_Q * rho * one_m_exp_alpha;
    psi_out[Mesh::HORZ] = Q + psi_in_V_minus_Q * one_m_exp_alpha;
  }
  else
  {
    double expf = exp_appx(-beta);
    double one_m_exp_beta = (1.0 - expf) / beta;
    psi_out[Mesh::VERT] = Q + psi_in_H_minus_Q * one_m_exp_beta;
    psi_out[Mesh::HORZ] = Q + psi_in_V_minus_Q * one_m_exp_beta / rho
                            + psi_in_H_minus_Q * (1.0 - 1.0/rho) * expf;
  }

  // Compute cell center flux.
  double psi_center = Q - (psi_out[Mesh::VERT] - psi_in[Mesh::VERT])/alpha
                        - (psi_out[Mesh::HORZ] - psi_in[Mesh::HORZ])/beta;

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi) psi[cell] = psi_center;

}

/// \todo It might be worth finding a faster exponential.
inline double Equation_SC_2D::exp_appx(double x)
{
//  double val = 1.0
//             + x
//             + 0.5 * x * x
//             + 0.1666666666666667 * x * x * x;
  double val = std::exp(x);
  return val;
}

} // end namespace detran

#endif /* detran_EQUATION_SC_2D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_2D.i.hh
//---------------------------------------------------------------------------//
