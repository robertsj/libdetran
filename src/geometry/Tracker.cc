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

  using std::cout;
  using std::endl;

  //-------------------------------------------------------------------------//
  // CREATE TRACK DATABASE AND INITIALIZE TRACKS
  //-------------------------------------------------------------------------//

  // Doing first two quadrants (eta > 0)
  d_trackdb = new TrackDB(2 * d_number_azimuths,
                          d_mesh->number_cells(),
                          d_quadrature);

  // Domain width (for x and y)
  double width = d_mesh->total_width_x();

  // Loop over all azimuths in the first two quadrants
  for (int a = 0; a < 2 * d_number_azimuths; a++)
  {

    // Add angular information.
    d_trackdb->setup_angle(a,
                           d_quadrature->cos_phi(a),
                           d_quadrature->sin_phi(a),
                           width*d_quadrature->spacing(a));

    for (int t = 0; t < d_quadrature->number_tracks(a); t++)
    {

      // Define the entrance and exit points on the actual mesh.
      // \todo This assumes a square mesh.
      Point enter = width * d_quadrature->enter(a, t);
      Point exit = width * d_quadrature->exit(a, t);

      // Create new track
      SP_track track(new Track(enter, exit));

      // Add the track to the database.
      d_trackdb->add_track(a, track);

    } // end track

  } // end angle


  // Do the actual track generation.
  generate_tracks();

}

void Tracker::normalize()
{
  vec_dbl volume(d_mesh->number_cells(), 0.0);
  for (int i = 0; i < volume.size(); i++)
    volume[i] = d_mesh->volume(i);
  // Normalize track lengths by volume.
  d_trackdb->normalize(volume);
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

  bool db = false;

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
    if (db) cout << "AZIMUTH = " << a << endl;
    if (db) cout << "    number tracks = " << d_trackdb->number_tracks_angle(a) << endl;

    // Loop through tracks
    for (int t = 0; t < d_trackdb->number_tracks_angle(a); t++)
    {
      if (db) cout << "    TRACK = " << t << endl;


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
      if (db) cout << "        tan_phi = " << tan_phi << endl;
      if (db) cout << "        sin_phi = " << sin_phi << endl;
      if (db) cout << "        cos_phi = " << cos_phi << endl;
      if (db) cout << "         length = " << track_length << endl;

      // Find the starting cell
      int IJ[] = {0, 0};
      find_starting_cell(enter, tan_phi, IJ);
      int I = IJ[0];
      int J = IJ[1];

      // Get segments
      double d_to_x = 0;
      double d_to_y = 0;
      p = enter;

      if (db) cout << "          ENTER = " << enter <<  endl;
      if (db) cout << "           EXIT = " << exit << endl;
      if (db) cout << "              I = " << I << endl;
      if (db) cout << "              J = " << J << endl;


      Assert(I <= d_mesh->number_cells_x());
      Assert(J <= d_mesh->number_cells_y());


      int count = 0;
      while (1)
      {

        if (db) cout << "        SEGMENT = " << count << endl;
        if (db) cout << "              I = " << I << endl;
        if (db) cout << "              J = " << J << endl;

        if (tan_phi > 0)
          d_to_x = d_x[I + 1] - p.x();
        else
          d_to_x = p.x() - d_x[I];

        d_to_y = d_y[J + 1] - p.y();

        // Flat source region.
        int region = d_mesh->index(I, J);
        if (db) cout << "            reg = " << region << endl;
        if (db) cout << "            d2x = " << d_to_x << endl;
        if (db) cout << "            d2y = " << d_to_y << endl;

        // Segment length
        double length = 0.0;

        double temp = std::abs(d_to_x * tan_phi) - d_to_y;

        if (db) cout << "           temp = " << temp << endl;
        if (std::abs(temp) > 1e-12 and temp > 0.0)
        {
          // I hit the top
          p = Point(p.x() + d_to_y / tan_phi, d_y[++J]);
          length = d_to_y / sin_phi;
          if (db) cout << "                NEW POINT 1 = " << p << endl;
        }
        else if (std::abs(temp) > 1e-12 and temp < 0.0)
        {
          // I hit the side
          if (tan_phi > 0.0)
            p = Point(d_x[++I], d_to_x * std::abs(tan_phi) + p.y());
          else
            p = Point(d_x[I--], d_to_x * std::abs(tan_phi) + p.y());
          length = d_to_x / std::abs(cos_phi);
          if (db) cout << "                NEW POINT 2 = " << p << endl;
        }
        else
        {
          // I cross through a corner
          if (tan_phi > 0.0)
            p = Point(d_x[++I], d_y[++J]);
          else
            p = Point(d_x[I--], d_y[++J]);
          length = d_to_y / sin_phi;
          if (db) cout << "                NEW POINT 3 = " << p << endl;
        }

        if (db) cout << "            len = " << length << endl;
        // Add a segment
        track->add_segment(Segment(region, length));

        // Check to see if we've left.
        if (I == -1 or I == d_x.size()-1 or J == d_y.size()-1)
        {
          double temp = distance(enter, p);

          if (db) cout << " lengths: " << temp <<  " " << track_length << " " << enter << " " << p << endl;
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

    if ( std::abs(d_x[i]-enter.x()) < 1e-10 )
    {
      if (tan_phi < 0.0) i--;
      break;
    }
    if (d_x[i] >= enter.x())
    {
      i--;
      break;
    }
  }
  if (i == d_x.size())
  {
    i--; i--;
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
  if (j == d_y.size())
  {
    j--; j--;
  }

  IJ[0] = i;
  IJ[1] = j;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file Tracker.cc
//---------------------------------------------------------------------------//
