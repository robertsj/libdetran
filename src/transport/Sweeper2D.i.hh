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

namespace detran
{

void Sweeper2D::sweep(moments_type &phi, moments_type &source)
{

  enum FACE2D
  {
    HORZ, VERT
  };

  // Reset the flux moments
  phi.assign(phi.size(), 0.0);

  // Temporary edge fluxes
  Equation::face_flux_2d psi_in  = {0.0, 0.0};
  Equation::face_flux_2d psi_out = {0.0, 0.0};

  // Sweep over all octants
  for (int o = 0; o < 4; o++)
  {

    d_equation->setup_octant(o);

    // Sweep over all angles
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {

      d_equation->setup_angle(a);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi)
      {
        psi = d_state->psi(o, a, d_g);
      }

      // Get boundary fluxes
      // psi_v
      // psi_h

      // Sweep over all y
      for (int j = 0; j < d_mesh->number_cells_y(); j++)
      {

        psi_out[VERT] = d_psi_v[j];

        for (int i = 0; i < d_mesh->number_cells_x(); i++)
        {

          psi_in[HORZ] = d_psi_h[i];
          psi_in[VERT] = psi_out[VERT];

          d_equation->solve(i, j, 0, source, psi_in, psi_out, phi, psi);

          // Save the horizontal flux.
          d_psi_h[i] = psi_out[HORZ];

        } // end x loop

        // Save the vertical flux.
        d_psi_v[j] = psi_out[VERT];

      } // end y loop

      // Update boundary
      // d_boundary->set???

    } // end angle loop

  } // end octant loop

  return;
}

}

#endif /* SWEEPER2D_I_HH_ */
