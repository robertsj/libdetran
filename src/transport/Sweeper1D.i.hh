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

#include "Equation_DD_1D.hh"

#include <iostream>

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper1D<EQ>::sweep(moments_type &phi)
{

  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

  // Temporary edge fluxes
  Equation<_1D>::face_flux_type psi_in = 0.0;
  Equation<_1D>::face_flux_type psi_out = 0.0;

  // Sweep over all octants
  for (int o = 0; o < 2; o++)
  {

    d_equation.setup_octant(o);

    // Sweep over all angles
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      SweepSource<_1D>::sweep_source_type source;
        source = d_sweepsource->source(d_g, o, a);

      // Setup equations for this angle.
      d_equation.setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(o, a, d_g);

      // Update (angle-wise)
      d_boundary->update(d_g, o, a);

      // Get boundary fluxes (are member v's needed?)
      psi_out =
        (*d_boundary)(d_face_index[o][Mesh::VERT][Boundary_T::IN], o, a, d_g);

      // Sweep over all cells.
      for (int ii = 0; ii < d_mesh->number_cells_x(); ii++)
      {
        // Get actual index.
        int i = index(o, 1, ii);

        psi_in = psi_out;

        d_equation.solve(i, 0, 0, source, psi_in, psi_out, phi, psi);

        // INSERT ACCELERATION MESH STUFF HERE
//        std::cout << " o = " << o << " a = " << a
//                  << " i = " << i << " psi_out = " << psi_out
//                  << " phi = " << phi[i] << std::endl;
      } // end x loop

      // Update boundary.
      (*d_boundary)(d_face_index[o][Mesh::VERT][Boundary_T::OUT], o, a, d_g) =
        psi_out;

    } // end angle loop

  } // end octant loop

  return;
}

// Instantiate
//template class Sweeper1D<Equation_SD_1D>;
template class Sweeper1D<Equation_DD_1D>;
//template class Sweeper1D<Equation_SC_1D>;

} // end namespace detran

#endif /* SWEEPER1D_I_HH_ */
