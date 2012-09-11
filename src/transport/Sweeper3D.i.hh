//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper3D.i.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper3D inline member definitions.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER3D_I_HH_
#define SWEEPER3D_I_HH_

#include "discretization/Equation_DD_3D.hh"
#include <iostream>
#ifdef DETRAN_ENABLE_OPENMP
#include <omp.h>
#endif

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper3D<EQ>::sweep(moments_type &phi)
{
  using std::cout;
  using std::endl;

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
  SweepSource<_3D>::sweep_source_type source(d_mesh->number_cells(), 0.0);

  // Sweep over all octants
  for (int o = 0; o < 8; o++)
  {

    equation.setup_octant(o);

    // Sweep over all angles
    #pragma omp for
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      d_sweepsource->source(d_g, o, a, source);

      // Setup equations for this angle.
      equation.setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(o, a, d_g);

      // Update the boundary for this angle.
      if (d_update_boundary) d_boundary->update(d_g, o, a);

      // Get boundary fluxes.
      boundary_flux_type psi_yz = (*d_boundary)
        (d_face_index[o][Mesh::YZ][Boundary_T::IN], o, a, d_g);
      boundary_flux_type psi_xz = (*d_boundary)
        (d_face_index[o][Mesh::XZ][Boundary_T::IN], o, a, d_g);
      boundary_flux_type psi_xy = (*d_boundary)
        (d_face_index[o][Mesh::XY][Boundary_T::IN], o, a, d_g);

      // Temporary edge fluxes.
      Equation<_3D>::face_flux_type psi_in  = { 0.0, 0.0, 0.0 };
      Equation<_3D>::face_flux_type psi_out = { 0.0, 0.0, 0.0 };

      // Sweep over all z
      for (int kk = 0; kk < d_mesh->number_cells_z(); kk++)
      {
        // Get actual index.
        int k = index(o, 3, kk);

        // Sweep over all y
        for (int jj = 0; jj < d_mesh->number_cells_y(); jj++)
        {
          // Get actual index.
          int j = index(o, 2, jj);

          psi_out[Mesh::YZ] = psi_yz[k][j];

          for (int ii = 0; ii < d_mesh->number_cells_x(); ii++)
          {
            // Get actual index.
            int i = index(o, 1, ii);

            psi_in[Mesh::YZ] = psi_out[Mesh::YZ];
            psi_in[Mesh::XZ] = psi_xz[k][i];
            psi_in[Mesh::XY] = psi_xy[j][i];

            // Solve.
            equation.solve(i, j, k, source, psi_in, psi_out, phi, psi);

            // Save the horizontal flux.
            psi_xz[k][i] = psi_out[Mesh::XZ];
            psi_xy[j][i] = psi_out[Mesh::XY];

            // INSERT ACCELERATION MESH STUFF HERE

          } // end x loop

          // Save the vertical flux.
          psi_yz[k][j] = psi_out[Mesh::YZ];

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
//template class Sweeper3D<Equation_SD_3D>;
template class Sweeper3D<Equation_DD_3D>;
//template class Sweeper3D<Equation_SC_3D>;

} // end namespace detran

#endif /* SWEEPER3D_I_HH_ */
