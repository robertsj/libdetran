//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper3D.t.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper 3D specialization.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER3D_T_HH_
#define SWEEPER3D_T_HH_

namespace detran
{

// Constructor setup.
template<>
  void
  Sweeper<_3D>::setup(SP_material material)
  {

    // Set up face/octant map.
    d_face_index.resize(8, vec2_int(3, vec_int(2, 0)));

    int inc[8][3] = {0, 2, 4,  1, 2, 4,  1, 3, 4,  0, 3, 4,
                     0, 2, 5,  1, 2, 5,  1, 3, 5,  0, 3, 5};
    int out[8][3] = {1, 3, 5,  0, 3, 5,  0, 2, 5,  1, 2, 5,
                     1, 3, 4,  0, 3, 4,  0, 2, 4,  1, 2, 4};
    for (int i = 0; i < 8; i++)
    {
      // octant     surface type    inout
      d_face_index[i][Mesh::YZ][Boundary_T::IN]  = inc[i][0];
      d_face_index[i][Mesh::XZ][Boundary_T::IN]  = inc[i][1];
      d_face_index[i][Mesh::XY][Boundary_T::IN]  = inc[i][2];
      d_face_index[i][Mesh::YZ][Boundary_T::OUT] = out[i][0];
      d_face_index[i][Mesh::XZ][Boundary_T::OUT] = out[i][1];
      d_face_index[i][Mesh::XY][Boundary_T::OUT] = out[i][2];
    }
  }

// Sweep.
template<>
  inline void
  Sweeper<_3D>::sweep(moments_type &phi)
  {

    // Reset the flux moments
    phi.assign(phi.size(), 0.0);

    // Temporary edge fluxes.
    Equation<_3D>::face_flux_type psi_in =
      { 0.0, 0.0, 0.0 };
    Equation<_3D>::face_flux_type psi_out =
      { 0.0, 0.0, 0.0 };

    // Temporary boundary fluxes.
    boundary_flux_type psi_yz;
    boundary_flux_type psi_xz;
    boundary_flux_type psi_xy;

    // Sweep over all octants
    for (int o = 0; o < 8; o++)
    {

      d_equation->setup_octant(o);

      // Sweep over all angles
      for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
      {

        // Get sweep source for this angle.
        SweepSource<_3D>::sweep_source_type source =
            d_sweepsource->source(d_g, o, a);

        // Setup equations for this angle.
        d_equation->setup_angle(a);

        // Get psi if update requested.
        State::angular_flux_type psi;
        if (d_update_psi) psi = d_state->psi(o, a, d_g);

        // Get boundary fluxes.
        psi_yz = (*d_boundary)(d_face_index[o][Mesh::YZ][Boundary_T::IN], o, a, d_g);
        psi_xz = (*d_boundary)(d_face_index[o][Mesh::XZ][Boundary_T::IN], o, a, d_g);
        psi_xy = (*d_boundary)(d_face_index[o][Mesh::XY][Boundary_T::IN], o, a, d_g);

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
              d_equation->solve(i, j, k, source, psi_in, psi_out, phi, psi);

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
        (*d_boundary)(d_face_index[o][Mesh::YZ][Boundary_T::OUT], o, a, d_g) = psi_yz;
        (*d_boundary)(d_face_index[o][Mesh::XZ][Boundary_T::OUT], o, a, d_g) = psi_xz;
        (*d_boundary)(d_face_index[o][Mesh::XY][Boundary_T::OUT], o, a, d_g) = psi_xy;

      } // end angle loop

    } // end octant loop

    return;
  }

// Instantiate
template class Sweeper<_3D> ;

} // end namespace detran

#endif /* SWEEPER3D_T_HH_ */
