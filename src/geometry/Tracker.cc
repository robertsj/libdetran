//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Tracker.cc
 *  @brief Tracker class member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Tracker.hh"
#include "angle/DPN.hh"
#include "utilities/SoftEquivalence.hh"
#include "utilities/MathUtilities.hh"
#include <iostream>

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Tracker::Tracker(SP_mesh mesh, SP_quadrature quadrature)
  : d_mesh(mesh)
  , d_quadrature(quadrature)
{
  Require(d_mesh);
  Insist(d_mesh->dimension() == 2, "Only 2D MOC is currently supported.");
  Require(d_quadrature);
  d_number_azimuths = d_quadrature->number_azimuths_octant();

  using std::cout;
  using std::endl;

  //--------------------------------------------------------------------------//
  // CREATE TRACK DATABASE AND INITIALIZE TRACKS
  //--------------------------------------------------------------------------//

  // Doing first two quadrants (eta > 0)
  d_trackdb = new TrackDB(2 * d_number_azimuths,
                          d_mesh->number_cells(),
                          d_quadrature);

  // Domain width and height (for x and y)
  double W = d_mesh->total_width_x();
  double H = d_mesh->total_width_y();

  double delta_max = 0.2; // 2 mm spacing
  using detran_utilities::vec_scale;
  // Loop over all azimuths in the first two quadrants
//  for (int a = 0; a < 2 * d_number_azimuths; a++)
//  {
//    double phi = d_quadrature->phi(a);
//    double cos_phi = d_quadrature->cos_phi(a);
//    double sin_phi = d_quadrature->sin_phi(a);
//    double tan_phi = sin_phi / cos_phi;
//    double cot_phi = 1.0 / tan_phi;
//    d_trackdb->setup_angle(a, cos_phi, sin_phi);
//
//    if (tan_phi < H / W)
//    {
//      // create the horizontal quadrature
//      n = number_points(phi, W, "H");
//      detran_angle::DPN q(n);
//      vec_dbl x  = q.cosines(0);
//      vec_dbl wt = q.weights();
//      vec_dbl X = x;
//      vec_scale(X, W);
//      vec_dbl wt_X = wt;
//      vec_scale(wt, W);
//      // create vertical that connects to horizontal
//      vec_dbl Y_h  = x;
//      vec_scale(Y_h, W * tan_phi);
//      vec_dbl wt_Y_h = wt;
//      vec_scale(wt_Y_h, W * tan_phi);
//      // create vertical that connects to vertical
//      n = number_points(phi, (H - W * tan_phi), "V");
//      detran_angle::DPN q2(n);
//      x  = q2.cosines(0);
//      wt = q2.weights();
//      vec_dbl Y_v = x;
//      vec_scale(Y_v, (H - W * tan_phi));
//      vec_dbl wt_Y_v = wt;
//      vec_scale(wt_Y_v, (H - W * tan_phi));
//      // piece together in order of incidence
//      for (size_t i = 0; i < X.size(); ++i)
//      {
//        Point r0(0.0,  Y_h[X.size()-i-1] + (H - W * tan_phi));
//        Point r1(X[i], H);
//        SP_track track(new Track(r0, r1, wt_X[i] * sin_phi));
//        d_trackdb->add_track(a, track);
//      }
//      for (size_t i = 0; i < Y_v.size(); ++i)
//      {
//        Point r0(0.0,  Y_v[i]);
//        Point r1(W,    Y_v[i] + W * tan_phi);
//        SP_track track(new Track(r0, r1, wt_Y_v[i] * cos_phi));
//        d_trackdb->add_track(a, track);
//      }
//      for (size_t i = 0; i < X.size(); ++i)
//      {
//        Point r0(X[i], 0.0);
//        Point r1(W,    Y_h[X.size()-i-1]);
//        SP_track track(new Track(r0, r1, wt_X[i] * sin_phi));
//        d_trackdb->add_track(a, track);
//      }
//    }
//    else
//    {
//      // create vertical quadrature
//      n = number_points(phi, delta_max, H, "V");
//      detran_angle::DPN q(n);
//      vec_dbl x  = q.cosines(0);
//      vec_dbl wt = q.weights();
//      vec_dbl Y = x;
//      vec_scale(Y, H);
//      vec_dbl wt_Y = wt;
//      vec_scale(wt, H);
//      // create horizontal that connects to vertical
//      vec_dbl X_v  = x;
//      vec_scale(X_v, H * cot_phi);
//      vec_dbl wt_X_v = wt;
//      vec_scale(wt_X_v, H * cot_phi);
//      // create horizontal that connects to horizontal
//      n = number_points(phi, delta_max, (W - H * cot_phi), "H");
//      detran_angle::DPN q2(n);
//      x  = q2.cosines(0);
//      wt = q2.weights();
//      vec_dbl X_h = x;
//      vec_scale(X_h, (H - H * cot_phi));
//      vec_dbl wt_X_h = wt;
//      vec_scale(wt_X_h, (H - H * cot_phi));
//      // piece together in order of incidence
//      for (size_t i = 0; i < X.size(); ++i)
//      {
//        Point r0(0.0,  Y_h[X.size()-i-1] + (H - W * tan_phi));
//        Point r1(X[i], H);
//        SP_track track(new Track(r0, r1, wt_X[i] * sin_phi));
//        d_trackdb->add_track(a, track);
//      }
//      for (size_t i = 0; i < Y_v.size(); ++i)
//      {
//        Point r0(0.0,  Y_v[i]);
//        Point r1(W,    Y_v[i] + W * tan_phi);
//        SP_track track(new Track(r0, r1, wt_Y_v[i] * cos_phi));
//        d_trackdb->add_track(a, track);
//      }
//      for (size_t i = 0; i < X.size(); ++i)
//      {
//        Point r0(X[i], 0.0);
//        Point r1(W,    Y_h[X.size()-i-1]);
//        SP_track track(new Track(r0, r1, wt_X[i] * sin_phi));
//        d_trackdb->add_track(a, track);
//      }
//    }
//
////    for (size_t t = 0; t < d_quadrature->number_tracks(a); t++)
////    {
////
////      // Define the entrance and exit points on the actual mesh.
////      // \todo This assumes a square mesh.
////      Point enter = width * d_quadrature->enter(a, t);
////      Point exit = width * d_quadrature->exit(a, t);
////
////      // Create new track
////      SP_track track(new Track(enter, exit, width*d_quadrature->spacing(a)));
////
////      // Add the track to the database.
////      d_trackdb->add_track(a, track);
////
////    } // end track
//
//  } // end angle

  // Do the actual track generation.
  generate_tracks();

}

/*
 * Basic model:
 *
 *          x
 *   .----------------.  H
 *   |      /         |
 *   |     .       -- |  y
 *   |    /     --    |
 *   |   .   --       |
 *   |  / -- __       |
 *   | --   /\        |
 *   |/       | phi   |
 *   .----------------.  0
 *   0                W
 *
 *   For line "--", note that
 *       y/W = tan(phi) --> y = W * tan(phi)
 *   and for line "-.-", note that
 *       H/x = tan(phi) --> x = H / tan(phi)
 */

//----------------------------------------------------------------------------//
void Tracker::generate_points(const double W,
                              const double H,
                              const size_t a)
{
//  using detran_utilities::vec_scale;
//  using detran_angle::DPN;
//  double phi     = d_quadrature->phi(a);
//  double cos_phi = d_quadrature->cos_phi(a);
//  bool octant = cos_phi < 0.0; // 0 or 1
//  //if (octant) cos_phi = std::abs(cos_phi);
//  double sin_phi = d_quadrature->sin_phi(a);
//  double tan_phi = sin_phi / cos_phi;
//  double cot_phi = 1.0 / tan_phi;
//
//  struct start_end
//  {
//    start_end(Point s = Point(), Point e = Point()) : start(s), end(e) {}
//    Point start;
//    Point end;
//  };
//  std::vector<start_end> points;
//
//  size_t dim[]   = {0, 1};
//  double width[] = {W, H};
//  double trig_phi = std::abs(tan_phi);
//  if (trig_phi < H / W)
//  {
//    dim[0]   = 1; dim[1] = 0;
//    width[0] = H; width[1] = W;
//    trig_phi = std::abs(cot_phi);
//  }
//  // create  quadrature for primary surface
//  DPN q(number_points(phi, width[0], dim[0]));
//  vec_dbl x = q.cosines(0);
//  vec_dbl w = q.weights();
//  vec_dbl P = x;
//  vec_scale(P, width[0]);
//  vec_dbl w_P = w;
//  vec_scale(w, width[0]);
//  // create secondary that connects to the primary
//  vec_dbl S_p  = x;
//  vec_scale(S_p, width[0] * trig_phi);
//  vec_dbl w_S_p = w;
//  vec_scale(w_S_p, width[0] * trig_phi);
//  // create secondary that connects to secondary
//  DPN q2(number_points(phi, (width[1] - width[0] * trig_phi), dim[1]));
//  x = q2.cosines(0);
//  w = q2.weights();
//  vec_dbl S_s = x;
//  vec_scale(S_s, (width[1] - width[0] * trig_phi));
//  vec_dbl w_S_s = w;
//  vec_scale(w_S_s, (width[1] - width[0] * trig_phi));
//  // add
//  for (size_t i = 0; i < P.size(); ++i)
//  {
//    double r[] = {0.0, 0.0};
//    r[dim[0]] =
//
//    Point r0(0.0,  S_p[P.size()-i-1] + (H - W * tan_phi));
//    Point r1(X[i], H);
//    SP_track track(new Track(r0, r1, wt_X[i] * sin_phi));
//    d_trackdb->add_track(a, track);
//  }
//  for (size_t i = 0; i < Y_v.size(); ++i)
//  {
//    Point r0(0.0,  Y_v[i]);
//    Point r1(W,    Y_v[i] + W * tan_phi);
//    SP_track track(new Track(r0, r1, wt_Y_v[i] * cos_phi));
//    d_trackdb->add_track(a, track);
//  }
//  for (size_t i = 0; i < X.size(); ++i)
//  {
//    Point r0(X[i], 0.0);
//    Point r1(W,    Y_h[X.size()-i-1]);
//    SP_track track(new Track(r0, r1, wt_X[i] * sin_phi));
//    d_trackdb->add_track(a, track);
//  }
}

//----------------------------------------------------------------------------//
Tracker::SP_tracker Tracker::Create(SP_mesh mesh, SP_quadrature quadrature)
{
  SP_tracker p(new Tracker(mesh, quadrature));
  return p;
}

//----------------------------------------------------------------------------//
void Tracker::normalize()
{
  vec_dbl volume(d_mesh->number_cells(), 0.0);
  for (size_t i = 0; i < volume.size(); i++)
    volume[i] = d_mesh->volume(i);
  // Normalize track lengths by volume.
  d_trackdb->normalize(volume);
}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

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
  for (size_t i = 0; i < d_x.size()-1; i++)
  {
    d_x[i + 1] = d_x[i] + d_mesh->dx(i);
  }
  for (size_t i = 0; i < d_y.size()-1; i++)
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
        if (std::abs(temp) > 1e-12 && temp > 0.0)
        {
          // I hit the top
          p = Point(p.x() + d_to_y / tan_phi, d_y[++J]);
          length = d_to_y / sin_phi;
          if (db) cout << "                NEW POINT 1 = " << p << endl;
        }
        else if (std::abs(temp) > 1e-12 && temp < 0.0)
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
        if (I == -1 || I == d_x.size()-1 || J == d_y.size()-1)
        {
          double temp = distance(enter, p);

          if (db) cout << " lengths: " << temp <<  " " << track_length << " " << enter << " " << p << endl;
          Ensure(detran_utilities::soft_equiv(temp, track_length));
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

//----------------------------------------------------------------------------//
Tracker::size_t Tracker::number_points(const double       phi,
                                       const double       L,
                                       const size_t       flag)
{
  Require(phi > 0.0);
  Require(phi < detran_utilities::pi * 0.5);
  double trig_phi = 0.0;
  if (!flag)
      trig_phi = std::sin(phi);
  else
      trig_phi = std::cos(phi);
  double tmp = std::min(d_delta_max, trig_phi * L);
  tmp = std::ceil((1.585/tmp)*(L*trig_phi) - 0.5851);
  size_t n = std::max((int)tmp, 3); // minimum of 3 tracks per surface
  return n;
}

//----------------------------------------------------------------------------//
void Tracker::segmentize(SP_track track)
{
  Point P0 = track->enter();
  Point P1 = track->exit();

  // Get all the x, y, and z intercepts
//  vec_int i = find(x, P0.x(), P0.x() > P1.x());
//  vec_int j = find(y, P0.y(), P0.x() > P1.x());
//  vec_int z = find(y, P0.y(), P0.x() > P1.x());
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Tracker.cc
//----------------------------------------------------------------------------//
