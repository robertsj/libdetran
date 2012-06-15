//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper1D.i.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper1D inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER1D_I_HH_
#define SWEEPER1D_I_HH_

// Detran
#include "Equation_DD_1D.hh"
#include "Equation_SD_1D.hh"

// System
#include <iostream>
#ifdef DETRAN_ENABLE_OPENMP
#include <omp.h>
#endif

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper1D<EQ>::sweep(moments_type &phi)
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
  SweepSource<_1D>::sweep_source_type source(d_mesh->number_cells(), 0.0);

  // Temporary edge fluxes
  Equation<_1D>::face_flux_type psi_in = 0.0;
  Equation<_1D>::face_flux_type psi_out = 0.0;

  // Sweep over all octants
  for (int o = 0; o < 2; o++)
  {

    // Setup equation for this octant.
    equation.setup_octant(o);

    // Sweep over all angles.
    #pragma omp for
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      d_sweepsource->source(d_g, o, a, source);

      // Setup equation for this angle.
      equation.setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(o, a, d_g);

      // Update the boundary for this angle.
      d_boundary->update(d_g, o, a);

      // Get boundary fluxes.
      psi_out =
        (*d_boundary)(d_face_index[o][Mesh::VERT][Boundary_T::IN], o, a, d_g);

      // Sweep over all cells.
      for (int ii = 0; ii < d_mesh->number_cells_x(); ii++)
      {
        // Get actual index.
        int i = index(o, 1, ii);

        // Set the incident cell surface flux.
        psi_in = psi_out;

        // Solve the equation in this cell.
        equation.solve(i, 0, 0, source, psi_in, psi_out, phi, psi);

        // ACCELERATION
        //if (d_acceleration) d_acceleration->tally(i, 0, 0, o, a, psi_out);

      } // end x loop

      // Update boundary.
      (*d_boundary)
        (d_face_index[o][Mesh::VERT][Boundary_T::OUT], o, a, d_g) = psi_out;

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
template class Sweeper1D<Equation_SD_1D>;
template class Sweeper1D<Equation_DD_1D>;
//template class Sweeper1D<Equation_SC_1D>;

} // end namespace detran

#endif /* SWEEPER1D_I_HH_ */
