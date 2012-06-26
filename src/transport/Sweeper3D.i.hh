//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper3D.i.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper3D inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER3D_I_HH_
#define SWEEPER3D_I_HH_

#include "Equation_DD_3D.hh"

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper3D<EQ>::sweep(moments_type &phi)
{

  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

  // Temporary edge fluxes.
  Equation<_3D>::face_flux_type psi_in  = { 0.0, 0.0, 0.0 };
  Equation<_3D>::face_flux_type psi_out = { 0.0, 0.0, 0.0 };

  // Temporary boundary fluxes.
  boundary_flux_type psi_yz;
  boundary_flux_type psi_xz;
  boundary_flux_type psi_xy;

  // Initialize equation and setup for this group.
  Equation_T equation(d_mesh, d_material, d_quadrature, d_update_psi);
  equation.setup_group(d_g);

  // Initialize discrete sweep source vector.
  SweepSource<_3D>::sweep_source_type source(d_mesh->number_cells(), 0.0);

  // Sweep over all octants
  for (int o = 0; o < 8; o++)
  {

    equation.setup_octant(o);

    // Sweep over all angles
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      d_sweepsource->source(d_g, o, a, source);

      // Setup equations for this angle.
      equation.setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(o, a, d_g);

      if (d_update_boundary) d_boundary->update(d_g, o, a);

      // Get boundary fluxes.
      psi_yz =
        (*d_boundary)(d_face_index[o][Mesh::YZ][Boundary_T::IN], o, a, d_g);
      psi_xz =
        (*d_boundary)(d_face_index[o][Mesh::XZ][Boundary_T::IN], o, a, d_g);
      psi_xy =
        (*d_boundary)(d_face_index[o][Mesh::XY][Boundary_T::IN], o, a, d_g);

      // Sweep over all z
      for (int k = 0; k < d_mesh->number_cells_z(); k++)
      {
        // Sweep over all y
        for (int j = 0; j < d_mesh->number_cells_y(); j++)
        {
          psi_out[Mesh::YZ] = psi_yz[k][j];

          for (int i = 0; i < d_mesh->number_cells_x(); i++)
          {
            psi_in[Mesh::YZ] = psi_out[Mesh::YZ];
            psi_in[Mesh::XZ] = psi_xz[k][i];
            psi_in[Mesh::XY] = psi_xy[j][i];

            // Solve.
            equation.solve(i, j, k, source, psi_in, psi_out, phi, psi);

            // Save the horizontal flux.
            psi_xz[i][k] = psi_out[Mesh::XZ];
            psi_xy[i][j] = psi_out[Mesh::XY];

            // INSERT ACCELERATION MESH STUFF HERE

          } // end x loop

          // Save the vertical flux.
          psi_yz[j][k] = psi_out[Mesh::YZ];

        } // end y loop
      } // end z loop

      // Update boundary
      (*d_boundary)
        (d_face_index[o][Mesh::YZ][Boundary_T::OUT], o, a, d_g) = psi_yz;
      (*d_boundary)
        (d_face_index[o][Mesh::XZ][Boundary_T::OUT], o, a, d_g) = psi_xz;
      (*d_boundary)
        (d_face_index[o][Mesh::XY][Boundary_T::OUT], o, a, d_g) = psi_xy;

      // Angular flux update
      if (d_update_psi) d_state->psi(d_g, o, a) = psi;

    } // end angle loop

  } // end octant loop

  #pragma omp master
  {
    d_number_sweeps++;
  }
  return;
}

// Instantiate
//template class Sweeper3D<Equation_SD_3D>;
template class Sweeper3D<Equation_DD_3D>;
//template class Sweeper3D<Equation_SC_3D>;

} // end namespace detran

#endif /* SWEEPER3D_I_HH_ */
