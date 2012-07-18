//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryMOC.cc
 * \brief  BoundaryMOC member definitions.
 * \author Jeremy Roberts
 * \date   Jun 26, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "BoundaryMOC.hh"
#include "VacuumMOC.hh"
#include "ReflectiveMOC.hh"

// System
#include <iostream>

namespace detran
{

// Constructor.
template<class D>
BoundaryMOC<D>::BoundaryMOC(SP_input        input,
                            SP_mesh         mesh,
                            SP_quadrature   quadrature)
  : Base(input, mesh)
  , d_quadrature(quadrature)
  , d_boundary_flux(d_number_groups,
                    vec3_dbl(quadrature->number_angles(),
                             vec2_dbl(2)))
  , d_bc(2*D::dimension)
{
  Require(d_quadrature);

  // Allocated the flux container.
  initialize();

  // Setup indexing.
  setup_indices();
  setup_side_indices();

  // Create boundary conditions.
  std::vector<std::string> names(6);
  names[Mesh::LEFT]   = "bc_left";
  names[Mesh::RIGHT]  = "bc_right";
  names[Mesh::BOTTOM] = "bc_bottom";
  names[Mesh::TOP]    = "bc_top";
  names[Mesh::SOUTH]  = "bc_south";
  names[Mesh::NORTH]  = "bc_north";

  // Assign boundary conditions.
  for(int side = 0; side < 2*D::dimension; side++)
  {
    // Vacuum is default.
    std::string type = "vacuum";
    if (input->check(names[side]))
    {
      type = input->get<std::string>(names[side]);
    }
    if (type == "vacuum")
    {
      d_bc[side] = new VacuumMOC<D>((*this), side, input, mesh, quadrature);
    }
    else if (type == "reflect")
    {
      d_bc[side] = new ReflectiveMOC<D>((*this), side, input, mesh, quadrature);
      d_is_reflective[side] = true;
      d_has_reflective = true;
    }
    else
    {
      type.append(" is not a supported bc type.");
      THROW(type);
      break;
    }
  }

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

template<class D>
void BoundaryMOC<D>::initialize()
{

  // The boundary flux array is already
  // sized to [groups, angles, 2], where 2 is for inout.

  // \todo Sloppy hack; need to specialize boundary so that the
  //       quadrature is specialized, too
  QuadratureMOC::SP_quadrature q;
  q = (QuadratureMOC::SP_quadrature) d_quadrature;

  for (int g = 0; g < d_number_groups; g++)
  {
    for (int o = 0; o < 4; o++)
    {
      for (int a = 0; a < q->number_angles_octant(); a++)
      {
        int azimuth = q->azimuth(a);
        int angle = q->index(o, a);
        d_boundary_flux[g][angle][IN].resize(q->number_tracks(azimuth), 0.0);
        d_boundary_flux[g][angle][OUT].resize(q->number_tracks(azimuth), 0.0);
      }
    }
  }

}

// \todo This is really ugly.  If time, clean up.
template<class D>
void BoundaryMOC<D>::setup_indices()
{
  using std::cout;
  using std::endl;

  // Each octant-angle-track is matched individually.  The indexing
  // is used as [angle, track]->(angle, track).

  int no = d_quadrature->number_octants();
  int na = d_quadrature->number_azimuths_octant();

  d_feed_into.resize(no, vec3_int(na));
  d_feed_from.resize(no, vec3_int(na));

  // Mirrored octants
  int oct[4][2] = {3, 1,
                   2, 0,
                   3, 1,
                   2, 0};

  // Loop over all angles
  for (int o = 0; o < no; o++)
  {
    for (int a = 0; a < na; a++)
    {

      // Number of exit intercepts on horizontal and vertical surface.
      int nx = d_quadrature->number_exit(a, 1);
      int ny = d_quadrature->number_exit(a, 0);

      // Size the arrays.
      d_feed_into[o][a].resize(d_quadrature->number_tracks(a), vec_int(3, 0));
      d_feed_from[o][a].resize(d_quadrature->number_tracks(a), vec_int(3, 0));

      // Assign the feed into index.
      int t_max = nx;
      if (o > 1) t_max = ny;
      for (int t = 0; t < t_max; t++)
      {
        // Example:  octant 0 into the top surface reflects into octant 3,
        //           same angle, but with an offset of ny (where ny is the number
        //           of tracks of the incoming angle that start on the y axis)
        int t2 = t;
        if (o > 1) t2 = t + nx;
        d_feed_into[o][a][t][0] = oct[o][0];
        d_feed_into[o][a][t][1] = a;
        d_feed_into[o][a][t][2] = t2;
      }
      t_max = ny;
      if (o > 1) t_max = nx;
      int offset = nx;
      if (o > 1) offset = ny;
      for (int t = 0; t < t_max; t++)
      {
        int t2 = t;
        if (o > 1) t2 = t + ny;
        d_feed_into[o][a][t + offset][0] = oct[o][1];
        d_feed_into[o][a][t + offset][1] = a;
        d_feed_into[o][a][t + offset][2] = t2;
      }

    }
  }

  // Loop over all angles
  for (int o = 0; o < no; o++)
  {
    for (int a = 0; a < na; a++)
    {

      // Number of exit intercepts on horizontal and vertical surface.
      int nx = d_quadrature->number_enter(a, 0);
      int ny = d_quadrature->number_enter(a, 1);

      // Assign the feed into index.
      for (int t = 0; t < nx + ny; t++)
      {
        int oo = d_feed_into[o][a][t][0];
        int aa = d_feed_into[o][a][t][1];
        int tt = d_feed_into[o][a][t][2];
        d_feed_from[oo][aa][tt][0] = o;
        d_feed_from[oo][aa][tt][1] = a;
        d_feed_from[oo][aa][tt][2] = t;
      }

    }
  }

//  // output side indices
//  for (int o = 0; o < d_feed_from.size(); o++)
//  {
//    cout << " OCTANT: " << o << endl;
//    for (int a = 0; a < d_feed_from[o].size(); a++)
//    {
//      cout << "   ANGLE: " << a << endl;
//      for (int t = 0; t < d_feed_from[o][a].size(); t++)
//      {
//        cout << "     TRACK: " << t << endl;
//        cout << "       FEEDS FROM " << endl;
//        cout << "           octant =  " << d_feed_from[o][a][t][0] << endl;
//        cout << "            angle =  " << d_feed_from[o][a][t][1] << endl;
//        cout << "            track =  " << d_feed_from[o][a][t][2] << endl;
//      }
//    }
//  }

//  // output side indices
//  for (int o = 0; o < d_feed_into.size(); o++)
//  {
//    cout << " OCTANT: " << o << endl;
//    for (int a = 0; a < d_feed_into[o].size(); a++)
//    {
//      cout << "   ANGLE: " << a << endl;
//      for (int t = 0; t < d_feed_into[o][a].size(); t++)
//      {
//        cout << "     TRACK: " << t << endl;
//        cout << "       FEEDS INTO " << endl;
//        cout << "           octant =  " << d_feed_into[o][a][t][0] << endl;
//        cout << "            angle =  " << d_feed_into[o][a][t][1] << endl;
//        cout << "            track =  " << d_feed_into[o][a][t][2] << endl;
//      }
//    }
//  }

}

// \todo I need refactoring!
template<class D>
void BoundaryMOC<D>::setup_side_indices()
{
  using std::cout;
  using std::endl;


  // Incident octants; left to right from incident perspective.
  int oct[4][2] = {0, 3,
                   2, 1,
                   1, 0,
                   3, 2};

  d_side_index.resize(4);

  vec_int triplet(3, 0);

  for (int side = 0; side < 2; side++)
  {

    // Octant 0.
    int o       = oct[side][0];
    int a_start = d_quadrature->number_azimuths_octant() - 1;
    int a_end   = 0;
    for (int a = a_start; a >= a_end; a--)
    {
      int azimuth = d_quadrature->azimuth(a);
      int ny      = d_quadrature->number_enter(azimuth, 1);
      for (int t = 0; t < ny; t++)
      {
        triplet[0] = o;
        triplet[1] = a;
        triplet[2] = t;
//        cout << " side = " << side << endl
//             << "      o = " << o
//             << "        a = " << a
//             << "          t = " << t << endl;
        d_side_index[side].push_back(triplet);
      }
    }

    // Octant 1.
    o       = oct[side][1];
    a_start = 0;
    a_end   = d_quadrature->number_azimuths_octant();
    for (int a = a_start; a < a_end; a++)
    {
      int azimuth = d_quadrature->azimuth(a);
      int ny = d_quadrature->number_enter(azimuth, 1);
      for (int t = ny - 1; t >= 0; t--)
      {
        triplet[0] = o;
        triplet[1] = a;
        triplet[2] = t;
        d_side_index[side].push_back(triplet);
      }
    }

  } // side loop 1

  for (int side = 2; side < 4; side++)
  {

    // Octant 0.
    int o       = oct[side][0];
    int a_start = 0;
    int a_end   = d_quadrature->number_azimuths_octant();
    for (int a = a_start; a < a_end; a++)
    {
      int azimuth = d_quadrature->azimuth(a);
      int nx = d_quadrature->number_enter(azimuth, 0);
      int ny = d_quadrature->number_enter(azimuth, 1);
      int n  = nx + ny;
      for (int t = n - 1; t >= ny; t--)
      {
        triplet[0] = o;
        triplet[1] = a;
        triplet[2] = t;
        d_side_index[side].push_back(triplet);
      }
    }

    // Octant 1.
    o       = oct[side][1];
    a_start = d_quadrature->number_azimuths_octant() - 1;
    a_end   = 0;
    for (int a = a_start; a >= a_end; a--)
    {
      int azimuth = d_quadrature->azimuth(a);
      int nx = d_quadrature->number_enter(azimuth, 0);
      int ny = d_quadrature->number_enter(azimuth, 1);
      int n  = nx + ny;
      for (int t = ny; t < n; t++)
      {
        triplet[0] = o;
        triplet[1] = a;
        triplet[2] = t;
        d_side_index[side].push_back(triplet);
      }
    }


  } // side loop 2

//  // output side indices
//  for (int s = 0; s < 4; s++)
//  {
//    cout << " >>> SIDE = " << s << endl;
//    for (int i = 0; i < d_side_index[s].size(); i++)
//    {
//      cout << "      o = " << d_side_index[s][i][0] << endl
//           << "        a = " << d_side_index[s][i][1] << endl
//           << "          t = " << d_side_index[s][i][2] << endl;
//    }
//  }

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.cc
//---------------------------------------------------------------------------//
