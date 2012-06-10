//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_SC_2D.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  Equation_SC_2D inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SC_2D_I_HH_
#define EQUATION_SC_2D_I_HH_

// System
#include <cmath>
#include <iostream>

namespace detran
{

inline void Equation_SC_2D::solve(int i,
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
  Require(i >= 0);
  Require(i < d_mesh->number_cells_x());
  Require(j >= 0);
  Require(j < d_mesh->number_cells_y());
  Require(k == 0);

  int cell = d_mesh->index(i, j);
  double sigma = d_material->sigma_t(d_mat_map[cell], d_g);
  double Q = source[cell] / sigma;
  double dx = d_mesh->dx(i);
  double dy = d_mesh->dy(j);
  double alpha = sigma * d_alpha[i];
  double beta = sigma * d_beta[j];
  double rho = alpha / beta; // \todo sigma cancels? precompute alpha/beta?

  // Compute outgoing face fluxes.
  if (rho <= 1.0)
  {
    double expf = exp_appx(-alpha);
    double one_m_exp_alpha = (1.0 - expf) / alpha;
    psi_out[Mesh::VERT] = Q + (psi_in[Mesh::VERT] - Q) * (1.0 - rho) * expf
                            + (psi_in[Mesh::HORZ] - Q) * rho * one_m_exp_alpha;
    psi_out[Mesh::HORZ] = Q + (psi_in[Mesh::VERT] - Q) * one_m_exp_alpha;
  }
  else
  {
    double expf = exp_appx(-beta);
    double one_m_exp_beta = (1.0 - expf) / beta;
    psi_out[Mesh::VERT] = Q + (psi_in[Mesh::HORZ] - Q) * one_m_exp_beta;
    psi_out[Mesh::HORZ] = Q + (psi_in[Mesh::VERT] - Q) * one_m_exp_beta / rho
                            + (psi_in[Mesh::HORZ] - Q) * (1.0 - 1.0/rho) * expf;
  }

  // Compute cell center flux.
  double psi_center = Q - (psi_out[Mesh::VERT] - psi_in[Mesh::VERT])/alpha
                        - (psi_out[Mesh::HORZ] - psi_in[Mesh::HORZ])/beta;

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi)
  {
    psi[cell] = psi_center;
  }
}

inline double Equation_SC_2D::exp_appx(double x)
{
  // exp(x) = exp(xo+d) ~ exp(xo)*(1 + d + d^2/2! + ... O(9))

  // Currently, there looks to be no benefit from the expansions.
  return std::exp(x);

  // Expressions from Maple CodeGeneration with optimization.
  double t1;
  double d;
  if (x <= 0.2)
  {
    t1 = 1.1051709180756476248; // exp(0.1)
    d  = x - 0.1;
  }
  else if (x <= 0.4)
  {
    t1 = 1.3498588075760031040; // exp(0.3)
    d  = x - 0.3;
  }
  else if (x <= 0.6)
  {
    t1 = 1.6487212707001281468; // exp(0.5)
    d  = x - 0.5;
  }
  else if (x <= 0.8)
  {
    t1 = 2.0137527074704765216; // exp(0.7)
    d  = x - 0.7;
  }
  else if (x <= 1.0)
  {
    t1 = 2.4596031111569496638; // exp(0.9)
    d  = x - 0.9;
  }
  else
  {
    // Otherwise, return the built in exp function.
    std::exp(x);
  }
  // 9th order
//  double t3 = d * d;
//  double t6 = t3 * d;
//  double t9 = t3 * t3;
//  double t21 = t9 * t9;
//  return t1 + t1 * d + t1 * t3 / 0.2e1 + t1 * t6 / 0.6e1 +
//         t1 * t9 / 0.24e2 + t1 * t9 * d / 0.120e3 +
//         t1 * t9 * t3 / 0.720e3 + t1 * t9 * t6 / 0.5040e4 +
//         t1 * t21 / 0.40320e5;
  // 8th order
//  double t3 = d * d;
//  double t6 = t3 * d;
//  double t9 = t3 * t3;
//  return t1 + t1 * d + t1 * t3 / 0.2e1 + t1 * t6 / 0.6e1 +
//         t1 * t9 / 0.24e2 + t1 * t9 * d / 0.120e3 +
//         t1 * t9 * t3 / 0.720e3 + t1 * t9 * t6 / 0.5040e4;
  // 7th order
  double t3 = d * d;
  double t9 = t3 * t3;
  return t1 + t1 * d + t1 * t3 / 0.2e1 +
         t1 * t3 * d / 0.6e1 + t1 * t9 / 0.24e2 +
         t1 * t9 * d / 0.120e3 + t1 * t9 * t3 / 0.720e3;


}

} // end namespace detran

#endif /* EQUATION_SC_2D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_2D.i.hh
//---------------------------------------------------------------------------//
