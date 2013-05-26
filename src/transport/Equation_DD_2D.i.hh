//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_DD_2D.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_DD_2D inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_DD_2D_I_HH_
#define detran_EQUATION_DD_2D_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline void Equation_DD_2D::solve(const size_t       i,
                                  const size_t       j,
                                  const size_t       k,
                                  moments_type      &source,
                                  face_flux_type    &psi_in,
                                  face_flux_type    &psi_out,
                                  moments_type      &phi,
                                  angular_flux_type &psi)
{
  // Preconditions.  (The client *must* set group and angles.)
  Require(k == 0);

  // Compute cell-center angular flux.
  int cell = d_mesh->index(i, j);
  double coef = 1.0 / (d_material->sigma_t(d_mat_map[cell], d_g) +
                       d_coef_x[i] + d_coef_y[j]);
  double psi_center = coef * (source[cell] +
                              d_coef_x[i] * psi_in[detran_geometry::Mesh::VERT] +
                              d_coef_y[j] * psi_in[detran_geometry::Mesh::HORZ] );

  // Compute outgoing fluxes.
  double two_psi_center = 2.0 * psi_center;
  psi_out[detran_geometry::Mesh::HORZ] =
    two_psi_center - psi_in[detran_geometry::Mesh::HORZ];
  psi_out[detran_geometry::Mesh::VERT] =
    two_psi_center - psi_in[detran_geometry::Mesh::VERT];

  // Compute flux moments.
  phi[cell] += d_quadrature->weight(d_angle) * psi_center;

  // Store angular flux if needed.
  if (d_update_psi) psi[cell] = psi_center;

}

} // end namespace detran

#endif /* detran_EQUATION_DD_2D_I_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_2D.i.hh
//---------------------------------------------------------------------------//
