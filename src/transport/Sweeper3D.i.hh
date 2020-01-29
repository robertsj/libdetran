//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper3D.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 24, 2012
 *  @brief  Sweeper3D inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef detran_SWEEPER3D_I_HH_
#define detran_SWEEPER3D_I_HH_

#include <iostream>
#ifdef DETRAN_ENABLE_OPENMP
#include <omp.h>
#endif

namespace detran
{

//---------------------------------------------------------------------------//
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

  // Reset the boundary flux tally
  if (d_tally) d_tally->reset(d_g);

  // Initialize discrete sweep source vector.
  SweepSource<_3D>::sweep_source_type source(d_mesh->number_cells(), 0.0);

  // Reference to boundary to simplify clutter.
  Boundary_T &b = *d_boundary;

  // Sweep over all octants
  for (size_t oo = 0; oo < 8; oo++)
  {
    size_t o = d_ordered_octants[oo];

    equation.setup_octant(o);

    // Sweep over all angles
    #pragma omp for
    for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
    {

      // Get sweep source for this angle.
      d_sweepsource->source(d_g, o, a, source);

      // Setup equations for this angle.
      equation.setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(d_g, o, a);

      // Update the boundary for this angle.
      if (d_update_boundary) b.update(d_g, o, a);

      // Get boundary fluxes.
      bf_type psi_yz = b(d_face_index[o][Mesh::YZ][Boundary_T::IN], o, a, d_g);
      bf_type psi_xz = b(d_face_index[o][Mesh::XZ][Boundary_T::IN], o, a, d_g);
      bf_type psi_xy = b(d_face_index[o][Mesh::XY][Boundary_T::IN], o, a, d_g);

      // Temporary edge fluxes.
      Equation<_3D>::face_flux_type psi_in  = { 0.0, 0.0, 0.0 };
      Equation<_3D>::face_flux_type psi_out = { 0.0, 0.0, 0.0 };

      // Sweep over all z
      int k  = d_space_ranges[o][2][0];
      int dk = d_space_ranges[o][2][1];
      for (size_t kk = 0; kk < d_mesh->number_cells_z(); ++kk, k += dk)
      {

        // Sweep over all y
        int j  = d_space_ranges[o][1][0];
        int dj = d_space_ranges[o][1][1];
        for (size_t jj = 0; jj < d_mesh->number_cells_y(); ++jj, j += dj)
        {

          psi_out[Mesh::YZ] = psi_yz[k][j];

          // Sweep over all x
          size_t i  = d_space_ranges[o][0][0];
          size_t di = d_space_ranges[o][0][1];
          for (size_t ii = 0; ii < d_mesh->number_cells_x(); ++ii, i += di)
          {
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
      b(d_face_index[o][Mesh::YZ][Boundary_T::OUT], o, a, d_g) = psi_yz;
      b(d_face_index[o][Mesh::XZ][Boundary_T::OUT], o, a, d_g) = psi_xz;
      b(d_face_index[o][Mesh::XY][Boundary_T::OUT], o, a, d_g) = psi_xy;

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

} // end namespace detran

#endif /* SWEEPER3D_I_HH_ */
