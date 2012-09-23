//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.i.hh
 * \author Jeremy Roberts
 * @date   Mar 24, 2012
 * \brief  Sweeper2D inline member definitions.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER2D_I_HH_
#define SWEEPER2D_I_HH_

#include "discretization/Equation_SD_2D.hh"
#include "discretization/Equation_DD_2D.hh"
#include "discretization/Equation_SC_2D.hh"
#include <iostream>
#ifdef DETRAN_ENABLE_OPENMP
#include <omp.h>
#endif

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper2D<EQ>::sweep(moments_type &phi)
{

  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

#ifdef DETRAN_ENABLE_OPENMP
  moments_type phi_local;
#else
  moments_type &phi_local = phi;
#endif

  #pragma omp parallel default(shared) private(phi_local)
  {

  // Initialize equation and setup for this group.
  Equation_T equation(d_mesh, d_material, d_quadrature, d_update_psi);
  equation.setup_group(d_g);

  // Reset the flux moments
  phi_local.resize(d_mesh->number_cells(), 0.0);

  // Initialize discrete sweep source vector.
  SweepSource<_2D>::sweep_source_type source(d_mesh->number_cells(), 0.0);

  // Sweep over all octants
  for (size_t o = 0; o < 4; o++)
  {

    // Setup equation for this octant.
    equation.setup_octant(o);

    // Get face indices
    const int face_V_i = d_face_index[o][Mesh::VERT][Boundary_T::IN];
    const int face_H_i = d_face_index[o][Mesh::HORZ][Boundary_T::IN];
    const int face_V_o = d_face_index[o][Mesh::VERT][Boundary_T::OUT];
    const int face_H_o = d_face_index[o][Mesh::HORZ][Boundary_T::OUT];

    // Sweep over all angles.
    #pragma omp for
    for (size_t a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      d_sweepsource->source(d_g, o, a, source);

      // Setup equation for this angle.
      equation.setup_angle(a);

      // Get psi if needed.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(d_g, o, a);

      // Update the boundary for this angle.
      if (d_update_boundary) d_boundary->update(d_g, o, a);

      // Get boundary fluxes.
      boundary_flux_type psi_v = (*d_boundary)(face_V_i, o, a, d_g);
      boundary_flux_type psi_h = (*d_boundary)(face_H_i, o, a, d_g);

      // Temporary edge fluxes.
      Equation<_2D>::face_flux_type psi_in  = {0.0, 0.0};
      Equation<_2D>::face_flux_type psi_out = {0.0, 0.0};

      // Sweep over all y.
      for (int jj = 0; jj < d_mesh->number_cells_y(); jj++)
      {
        // Get actual index.
        int j = index(o, 2, jj);

        // Note the index: incident boundaries are access from
        // "left to right" w/r to self.
        psi_out[Mesh::VERT] = psi_v[j];

        // Sweep over all x.
        for (int ii = 0; ii < d_mesh->number_cells_x(); ii++)
        {
          // Get actual index.
          int i = index(o, 1, ii);

          // Set the incident cell surface fluxes.
          psi_in[Mesh::HORZ] = psi_h[i];
          psi_in[Mesh::VERT] = psi_out[Mesh::VERT];

          // Solve the equation in this cell.
          equation.solve(i, j, 0, source, psi_in, psi_out, phi_local, psi);

          // Save the horizontal flux.
          psi_h[i] = psi_out[Mesh::HORZ];

          // INSERT ACCELERATION MESH STUFF HERE

        } // end x loop

        // Save the vertical flux.
        psi_v[j] = psi_out[Mesh::VERT];

      } // end y loop

      // Update boundary
      (*d_boundary)(face_V_o, o, a, d_g) = psi_v;
      (*d_boundary)(face_H_o, o, a, d_g) = psi_h;

      // Update the angular flux.
      if (d_update_psi) d_state->psi(d_g, o, a) = psi;

    } // end angle loop
    // end omp do

  } // end octant loop

#ifdef DETRAN_ENABLE_OPENMP
  // Sum local thread fluxes.
  #pragma omp critical
  {
    for (int i = 0; i < d_mesh->number_cells(); i++)
    {
      phi[i] += phi_local[i];
    }
  }
#endif

  } // end omp parallel

  #pragma omp master
  {
    d_number_sweeps++;
  }
  return;
}

// Instantiate
template class Sweeper2D<Equation_DD_2D>;
template class Sweeper2D<Equation_SC_2D>;
template class Sweeper2D<Equation_SD_2D>;

} // end namespace detran

#endif /* SWEEPER2D_I_HH_ */
