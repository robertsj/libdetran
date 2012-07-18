//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2DMOC.i.hh
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Sweeper2DMOC inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER2DMOC_I_HH_
#define SWEEPER2DMOC_I_HH_

// Detran
#include "Equation_SC_MOC.hh"

// System
#include <iostream>
#ifdef DETRAN_ENABLE_OPENMP
#include <omp.h>
#endif

namespace detran
{

// Sweep.
template <class EQ>
inline void Sweeper2DMOC<EQ>::sweep(moments_type &phi)
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
  SweepSource<_2D>::sweep_source_type source(d_mesh->number_cells(), 0.0);

  double psi_in  = 0;
  double psi_out = 0;

  QuadratureMOC::SP_quadrature q = d_quadrature;

  // Sweep over all octants.
  for (int o = 0; o < 4; o++)
  {
    //cout << "OCTANT: " << o << endl;

    // Setup equation for this octant.
    equation.setup_octant(o);

    // Sweep over all angles.
    #pragma omp for
    for (int a = 0; a < q->number_angles_octant(); a++)
    {
      //cout << "  ANGLE: " << a << endl;

      // Get the azimuth and polar indices within the octant.
      int azimuth = q->azimuth(a);
      int polar   = q->polar(a);
      //cout << "  AZIMUTH: " << azimuth << endl;
      //cout << "  POLAR: " << polar << endl;

      // Switch this back to one angle setup.
      equation.setup_azimuth(azimuth);
      equation.setup_polar(polar);

      // Switch the azimuth index to correct one for track access.
      bool track_reverse = false;
      switch (o)
      {
        case 0 :
          // do nothing
          break;
        case 1:
          azimuth += q->number_azimuths_octant();
          break;
        case 2:
          track_reverse = true;
          break;
        case 3:
          azimuth += q->number_azimuths_octant();
          track_reverse = true;
          break;
        default:
          THROW("WTF!!!");
          break;
      }

      // Get sweep source for this angle.
      d_sweepsource->source(d_g, o, a, source);

      // Get psi if update requested.
      State::angular_flux_type psi;
      if (d_update_psi) psi = d_state->psi(d_g, o, a);

      // Update the boundary for this angle.
      if (d_update_boundary) d_boundary->update(d_g, o, a);

      // Sweep over all tracks.
      for (int t = 0; t < d_tracks->number_tracks_angle(azimuth); t++)
      {
//        if (!track_reverse)
//          cout << "    TRACK: " << t << endl;
//        else
//          cout << "    TRACK: " << t << " in REVERSE" << endl;
        // Get track.
        SP_track track = d_tracks->track(azimuth, t);

        // *** LOAD THE BOUNDARY FLUX.
//        if (o == 0)
//          psi_out = 1.0 + 4*o + t;
//        else
        psi_out = (*d_boundary)(d_g, o, a, BoundaryMOC<_2D>::IN, t);

//        cout << " OCTANT = " << o << endl;
//        cout << "   ANGLE = " << a << endl;
//        cout << "     TRACK = " << t << endl;
//        cout << "       psi_in = " << psi_out << endl;

        // SN access
        // boundary_flux_type psi_v = (*d_boundary)
        //   (d_face_index[o][Mesh::VERT][Boundary_T::IN], o, a, d_g);
        // MOC access
        // double psi_v = (*d_boundary)(d_g, o, a, t);
        // --> the side doesn't matter for this access
        // --> ergo, the side really needn't be part of the storage
        // --> create index maps for which track is on a side, etc.


        //cout << "      PSI_IN: " << psi_out << endl;

        // Sweep all segments on the track.
        int s = 0;
        for (int ss = 0; ss < track->number_segments(); ss++)
        {
          s = ss;
          if (track_reverse) s = track->number_segments() - ss - 1;

          //cout << "      SEGMENT: " << s << endl;

          // Update track angular flux
          psi_in = psi_out;

          // Get segment region.
          int region = track->segment(s).region();

          // Get segment length.
          double length = track->segment(s).length();

          // Solve.
          equation.solve(region, length, source, psi_in, psi_out, phi_local, psi);

          //cout << "        PSI_OUT: " << psi_out << endl;

        } // end segment

        // *** UPDATE THE BOUNDARY WITH psi_out
        (*d_boundary)(d_g, o, a, BoundaryMOC<_2D>::OUT, t) = psi_out;

        //cout << "      psi_out = " << psi_out << endl;

      } // end track



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
template class Sweeper2DMOC<Equation_SC_MOC>;

} // end namespace detran

#endif /* SWEEPER2DMOC_I_HH_ */
