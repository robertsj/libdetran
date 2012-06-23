//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Tracker.cc
 * \brief  Tracker 
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "Tracker.hh"

// Utilities
#include "SoftEquivalence.hh"

// System
#include <iostream>

namespace detran
{

Tracker::Tracker(SP_mesh mesh, SP_quadrature quadrature)
  : d_mesh(mesh)
  , d_quadrature(quadrature)
{
  Require(d_mesh);
  Insist(d_mesh->dimension() == 2, "Only 2D MOC is currently supported.");
  Insist(soft_equiv(d_mesh->total_width_x(), d_mesh->total_width_y()),
    "MOC is currently supported only on square domains");
  Require(d_quadrature);
  d_number_azimuths = d_quadrature->number_azimuths_octant();

  //-------------------------------------------------------------------------//
  // CREATE TRACK DATABASE AND INITIALIZE TRACKS
  //-------------------------------------------------------------------------//

  // Doing first two quadrants (eta > 0)
  d_trackdb = new TrackDB(2 * d_number_azimuths);

  // Domain width (for x and y)
  double width = d_mesh->total_width_x();

  // Loop over all azimuths in the first two quadrants
  for (int a = 0; a < 2 * d_number_azimuths; a++)
  {

    for (int t = 0; t < d_quadrature->number_tracks(a); t++)
    {

      // Define the entrance and exit points on the actual mesh.
      // \todo This assumes a square mesh.
      Point enter = width * d_quadrature->enter(a, t);
      Point exit = width * d_quadrature->exit(a, t);

      std::cout << " a = " << a << " t = " << t
                << " enter = " << enter << " exit = " << exit << std::endl;

      // Create new track
      SP_track track(new Track(enter, exit));

      // Add the track to the database.
      d_trackdb->add_track(a, track);

    } // end track

  } // end angle


  // Do the actual track generation.
  generate_tracks();

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

void Tracker::generate_tracks()
{

  using std::cout;
  using std::endl;
  /*
   *  for all azimuth
   *    for all origins
   *      get track
   *      while not at exit
   *        ray trace the grid
   */

  // Create the mesh grid.
  d_x.resize(d_mesh->number_cells_x() + 1, 0.0);
  d_y.resize(d_mesh->number_cells_y() + 1, 0.0);
  for (int i = 0; i < d_x.size()-1; i++)
  {
    d_x[i + 1] = d_x[i] + d_mesh->dx(i);
  }
  for (int i = 0; i < d_y.size()-1; i++)
  {
    d_y[i + 1] = d_y[i] + d_mesh->dy(i);
  }

  // Loop through azimuths
  for (int a = 0; a < d_number_azimuths * 2; a++)
  {
    cout << "AZIMUTH = " << a << endl;
    cout << "    number tracks = " << d_trackdb->number_tracks_angle(a) << endl;

    // Loop through tracks
    for (int t = 0; t < d_trackdb->number_tracks_angle(a); t++)
    {
      cout << "    TRACK = " << t << endl;


      // Get the track
      SP_track track = d_trackdb->track(a, t);

      // Compute the track length
      Point enter = track->enter();
      Point exit  = track->exit();
      double track_length = distance(enter, exit);

      // Compute tangent of angle with respect to x
      Point p = exit - enter;
      double tan_phi = p.y() / p.x();
      double sin_phi = (exit.y()-enter.y()) / track_length;
      double cos_phi = (exit.x()-enter.x()) / track_length;
      cout << "        tan_phi = " << tan_phi << endl;
      cout << "        sin_phi = " << sin_phi << endl;
      cout << "        cos_phi = " << cos_phi << endl;
      cout << "         length = " << track_length << endl;

      // Find the starting cell
      int IJ[] = {0, 0};
      find_starting_cell(enter, tan_phi, IJ);
      int I = IJ[0];
      int J = IJ[1];

      Assert(I <= d_mesh->number_cells_x());
      Assert(J <= d_mesh->number_cells_y());

      // Get segments
      double d_to_x = 0;
      double d_to_y = 0;
      p = enter;

      cout << "          ENTER = " << enter <<  endl;
      cout << "           EXIT = " << exit << endl;
      cout << "              I = " << I << endl;
      cout << "              J = " << J << endl;
      int count = 0;
      while (1)
      {

        cout << "        SEGMENT = " << count << endl;
        cout << "              I = " << I << endl;
        cout << "              J = " << J << endl;

        if (tan_phi > 0)
          d_to_x = d_x[I + 1] - p.x();
        else
          d_to_x = p.x() - d_x[I];

        d_to_y = d_y[J + 1] - p.y();

        // Flat source region.
        int region = d_mesh->index(I, J);
        cout << "            reg = " << region << endl;
        cout << "            d2x = " << d_to_x << endl;
        cout << "            d2y = " << d_to_y << endl;

        // Segment length
        double length = 0.0;

        double temp = std::abs(d_to_x * tan_phi) - d_to_y;

        cout << "           temp = " << temp << endl;
        if (std::abs(temp) > 1e-10 and temp > 0.0)
        {
          // I hit the top
          p = Point(p.x() + d_to_y / tan_phi, d_y[++J]);
          length = d_to_y / sin_phi;
          cout << "                NEW POINT 1 = " << p << endl;
        }
        else if (std::abs(temp) > 1e-10 and temp < 0.0)
        {
          // I hit the side
          if (tan_phi > 0.0)
            p = Point(d_x[++I], d_to_x * std::abs(tan_phi) + d_y[J]);
          else
            p = Point(d_x[I--], d_to_x * std::abs(tan_phi) + d_y[J]);
          length = d_to_x / std::abs(cos_phi);
          cout << "                NEW POINT 2 = " << p << endl;
        }
        else
        {
          // I cross through a corner
          if (tan_phi > 0.0)
            p = Point(d_x[++I], d_y[++J]);
          else
            p = Point(d_x[I--], d_y[++J]);
          length = d_to_y / sin_phi;
          cout << "                NEW POINT 3 = " << p << endl;
        }

        cout << "            len = " << length << endl;
        // Add a segment
        track->add_segment(Segment(region, length));

        // Check to see if we've left.
        if (I == -1 or I == d_x.size()-1 or J == d_y.size()-1)
        {
          double temp = distance(enter, p);

          cout << " lengths: " << temp <<  " " << track_length << " " << enter << " " << p << endl;
          Ensure(soft_equiv(temp, track_length));
          break;
        }
        count++;

      } // segment loop

    } // track loop

  } // azimuth loop


}

// Find the indices of the cell we are to enter
void Tracker::find_starting_cell(Point enter, double tan_phi, int *IJ)
{
  int i = 0;
  int j = 0;

  // Going right.

  for (i = 1; i < d_x.size(); i++)
  {
    if (d_x[i] >= enter.x())
    {
      i--;
      break;
    }
  }
  //if (tan_phi < 0) i--;
  Assert(i >= 0);


  for (j = 1; i < d_y.size() - 1; j++)
  {
    if (d_y[j] > enter.y())
    {
      j--;
      break;
    }
  }
  IJ[0] = i;
  IJ[1] = j;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Tracker.cc
//---------------------------------------------------------------------------//
