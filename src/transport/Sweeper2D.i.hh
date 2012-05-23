//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.i.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper2D inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER2D_I_HH_
#define SWEEPER2D_I_HH_

#include "Equation_DD_2D.hh"
#include "Equation_SC_2D.hh"

#include <iostream>

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper2D<EQ>::sweep(moments_type &phi)
{

  d_equation.setup_group(d_g);

  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

  // Sweep over all octants
  for (int o = 0; o < 4; o++)
  {

    d_equation.setup_octant(o);

    // Sweep over all angles
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      SweepSource<_2D>::sweep_source_type source =
          d_sweepsource->source(d_g, o, a);

      // Setup equations for this angle.
      d_equation.setup_angle(a);

      // Get psi if needed.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(d_g, o, a);

      // Update (angle-wise)
      d_boundary->update(d_g, o, a);

      // Get boundary fluxes (are member v's needed?) \todo
      boundary_flux_type psi_v = (*d_boundary)
          (d_face_index[o][Mesh::VERT][Boundary_T::IN], o, a, d_g);
      boundary_flux_type psi_h = (*d_boundary)
          (d_face_index[o][Mesh::HORZ][Boundary_T::IN], o, a, d_g);

      // Temporary edge fluxes
      Equation<_2D>::face_flux_type psi_in  = {0.0, 0.0};
      Equation<_2D>::face_flux_type psi_out = {0.0, 0.0};

      // Sweep over all y
      for (int jj = 0; jj < d_mesh->number_cells_y(); jj++)
      {
        // Get actual index.
        int j = index(o, 2, jj);

        psi_out[Mesh::VERT] = psi_v[j];

        for (int ii = 0; ii < d_mesh->number_cells_x(); ii++)
        {
          // Get actual index.
          int i = index(o, 1, ii);

          psi_in[Mesh::HORZ] = psi_h[i];
          psi_in[Mesh::VERT] = psi_out[Mesh::VERT];

          // Solve the equation in this cell.
          d_equation.solve(i, j, 0, source, psi_in, psi_out, phi, psi);

          // Save the horizontal flux.
          psi_h[i] = psi_out[Mesh::HORZ];

          // INSERT ACCELERATION MESH STUFF HERE

        } // end x loop

        // Save the vertical flux.
        psi_v[j] = psi_out[Mesh::VERT];

      } // end y loop

      // Update boundary
      (*d_boundary)
          (d_face_index[o][Mesh::VERT][Boundary_T::OUT], o, a, d_g) = psi_v;
      (*d_boundary)
          (d_face_index[o][Mesh::HORZ][Boundary_T::OUT], o, a, d_g) = psi_h;

      // Angular flux update
      if (d_update_psi) d_state->psi(d_g, o, a) = psi;

    } // end angle loop

  } // end octant loop

  //STOP_PROFILER();

  return;
}

// Instantiate
//template class Sweeper1D<Equation_SD_2D>;
template class Sweeper2D<Equation_DD_2D>;
template class Sweeper2D<Equation_SC_2D>;

} // end namespace detran

#endif /* SWEEPER2D_I_HH_ */
