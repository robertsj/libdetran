//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper1D.t.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper 1D specialization.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER1D_T_HH_
#define SWEEPER1D_T_HH_

namespace detran
{

// Constructor setup.
template <>
inline void Sweeper<_1D>::setup(SP_material material)
{
  // Set up face/octant map.
  d_face_index.resize(2, vec2_int(1, vec_int(2, 0)));
  // octant / surface type / inout
  d_face_index[0][Mesh::VERT][Boundary_T::IN] = Mesh::LEFT;
  d_face_index[1][Mesh::VERT][Boundary_T::IN] = Mesh::RIGHT;
  d_face_index[0][Mesh::VERT][Boundary_T::OUT] = Mesh::RIGHT;
  d_face_index[1][Mesh::VERT][Boundary_T::OUT] = Mesh::LEFT;
}

// Sweep.
template <>
inline void Sweeper<_1D>::sweep(moments_type &phi)
{

  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

  // Temporary edge fluxes
  Equation<_1D>::face_flux_type psi_in = 0.0;
  Equation<_1D>::face_flux_type psi_out = 0.0;

  // Sweep over all octants
  for (int o = 0; o < 2; o++)
  {

    d_equation->setup_octant(o);

    // Sweep over all angles
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      SweepSource<_1D>::sweep_source_type source = d_sweepsource->source(d_g,
                                                                         o, a);

      // Setup equations for this angle.
      d_equation->setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi)
      {
        psi = d_state->psi(o, a, d_g);
      }

      // Get boundary fluxes (are member v's needed?)
      psi_out = (*d_boundary)(d_face_index[o][Mesh::VERT][Boundary_T::IN], o,
                              a, d_g);

      // Sweep over all cells.
      for (int i = 0; i < d_mesh->number_cells_x(); i++)
      {

        psi_in = psi_out;

        d_equation->solve(i, 0, 0, source, psi_in, psi_out, phi, psi);

        // INSERT ACCELERATION MESH STUFF HERE

      } // end x loop

      // Update boundary.
      (*d_boundary)(d_face_index[o][Mesh::VERT][Boundary_T::OUT], o, a, d_g) =
          psi_out;

    } // end angle loop

  } // end octant loop

  return;
}

// Instantiate
template class Sweeper<_1D> ;

} // end namespace detran

#endif /* SWEEPER1D_T_HH_ */
