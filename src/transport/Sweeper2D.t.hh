//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.t.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper 2D specialization.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER2D_T_HH_
#define SWEEPER2D_T_HH_

// Detran
#include "Equation_DD_2D.hh"

// System
#include <iostream>

namespace detran
{
using std::cout;
using std::endl;

template<>
void Sweeper<_2D>::setup(SP_material material)
{
  // Check to see if input has equation type. If not, default.
  std::string equation;
  if (d_input->check("equation"))
  {
    cout << "equation found." << endl;
    equation = d_input->get<std::string>("equation");
  }
  else
  {
    equation = "dd";
  }
  // Set the equation.
  if (equation == "dd")
  {
    d_equation =
        new Equation_DD_2D(d_mesh, material, d_quadrature, d_update_psi);
  }
  else
  {
    std::cout << " equation = " << equation << std::endl;
    THROW("Unsupported equation type");
  }

  // Set up face/octant map.
  d_face_index.resize(4, vec2_int(2, vec_int(2, 0)));
  // octant / surface type / inout
  for (int i = 0; i < 4; i++)
  {
    d_face_index[i][Mesh::VERT][Boundary_T::IN]  = i % 2 + Mesh::LEFT;
    d_face_index[i][Mesh::HORZ][Boundary_T::IN]  = i % 2 + Mesh::BOTTOM;
    d_face_index[i][Mesh::VERT][Boundary_T::OUT] = Mesh::RIGHT - (i + 1) % 2;
    d_face_index[i][Mesh::HORZ][Boundary_T::OUT] = Mesh::TOP   - (i + 1) % 2;
  }
}

template<>
inline void Sweeper<_2D>::sweep(moments_type &phi)
{
  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

  // Temporary edge fluxes
  Equation<_2D>::face_flux_type psi_in  = {0.0, 0.0};
  Equation<_2D>::face_flux_type psi_out = {0.0, 0.0};

  // Reference boundary fluxes for notational clarity.
  boundary_flux_type psi_v;
  boundary_flux_type psi_h;

  // Mesh loop bounds.  This is temporary: a cleaner approach will include adjoint.
  int ib[4][3] = { {0, d_mesh->number_cells_x(),  1},
                   {d_mesh->number_cells_x(), 0, -1},
                   {d_mesh->number_cells_x(), 0, -1},
                   {0, d_mesh->number_cells_x(),  1}};
  int jb[4][3] = { {0, d_mesh->number_cells_y(),  1},
                   {0, d_mesh->number_cells_y(),  1},
                   {d_mesh->number_cells_y(), 0, -1},
                   {d_mesh->number_cells_y(), 0, -1}};

  // Sweep over all octants
  for (int o = 0; o < 4; o++)
  {

    d_equation->setup_octant(o);

    // Sweep over all angles
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      // Get sweep source for this angle.
      SweepSource<_2D>::sweep_source_type source =
          d_sweepsource->source(d_g, o, a);

      // Setup equations for this angle.
      d_equation->setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi)
      {
        THROW(" why are you here? ");
        psi = d_state->psi(o, a, d_g);
      }

      // Get boundary fluxes (are member v's needed?)
      psi_v = (*d_boundary)
          (d_face_index[o][Mesh::VERT][Boundary_T::IN], o, a, d_g);
      psi_h = (*d_boundary)
          (d_face_index[o][Mesh::HORZ][Boundary_T::IN], o, a, d_g);

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

          d_equation->solve(i, j, 0, source, psi_in, psi_out, phi, psi);

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

    } // end angle loop

  } // end octant loop

  return;
}

// Instantiate
template class Sweeper<_2D>;

} // end namespace detran

#endif /* SWEEPER2D_T_HH_ */
