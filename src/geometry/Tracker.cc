//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Tracker.cc
 *  @brief Tracker class member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Tracker.hh"
#include "angle/QuadratureFactory.hh"
#include "utilities/SoftEquivalence.hh"
#include "utilities/MathUtilities.hh"
#include <iostream>

#define COUT(c) std::cout << c << std::endl;

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Tracker::Tracker(SP_db db, SP_quadrature q)
  : d_db(db)
  , d_quadrature(q)
  , d_X(0)
  , d_Y(0)
  , d_Z(0)
  , d_maximum_spacing(0.1)
  , d_spatial_quad_type("uniform")
{
  Require(d_db);
  Require(d_quadrature);

  if (d_db->check("tracker_maximum_spacing"))
  {
    d_maximum_spacing = d_db->get<double>("tracker_maximum_spacing");
    Assert(d_maximum_spacing > 0.0);
  }
  if (d_db->check("tracker_spatial_quad_type"))
  {
    d_spatial_quad_type = d_db->get<std::string>("tracker_spatial_quad_type");
  }

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
Tracker::SP_tracker Tracker::Create(SP_db db, SP_quadrature quadrature)
{
  SP_tracker p(new Tracker(db, quadrature));
  return p;
}

//----------------------------------------------------------------------------//
void Tracker::normalize()
{
  vec_dbl volume(d_mesh->number_cells(), 0.0);
  for (size_t i = 0; i < volume.size(); i++)
    volume[i] = d_mesh->volume(i);
  // Normalize track lengths by volume.
  d_tracks->normalize(volume);
}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

// 2-d mesh
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
  for (int a = 0; a < d_quadrature->number_azimuths_octant() * 2; a++)
  {
    if (db) cout << "AZIMUTH = " << a << endl;
    if (db) cout << "    number tracks = " << d_tracks->number_tracks_angle(a) << endl;

    // Loop through tracks
    for (int t = 0; t < d_tracks->number_tracks_angle(a); t++)
    {
      if (db) cout << "    TRACK = " << t << endl;


      // Get the track
      SP_track track = d_tracks->track(a, t);

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
Tracker::size_t Tracker::number_points(const double  phi,
                                       const double  L,
                                       const size_t  flag)
{
  Require(phi > 0.0);
  Require(phi < detran_utilities::pi * 0.5);
  double trig_phi = 0.0;
  if (!flag)
      trig_phi = std::sin(phi);
  else
      trig_phi = std::cos(phi);


  double tmp = std::min(d_maximum_spacing, trig_phi * L);


  COUT("phi " << phi)
  COUT("L " << L)
  COUT("flat " << flag)
  COUT("trig_phi " << trig_phi)
  COUT("max_d " << tmp)

  if (d_spatial_quad_type == "gl")
  {
    tmp = std::ceil((1.585/tmp)*(L*trig_phi) - 0.5851);
  }
  else
  {
    tmp = std::ceil(L * trig_phi / tmp);
  }

  COUT("n "  << tmp)


  size_t n = std::max((int)tmp, 3); // minimum of 3 tracks per surface
  COUT("computed = " << tmp << " actual = " << n)
  return n;

//  if (flag == "H") :
//      trig_phi = np.sin(phi)
//  else :
//      trig_phi = np.cos(phi)
//  max_d = min(max_d, trig_phi * L)
//  return max(int(np.ceil((1.585/max_d)*(L*trig_phi) - 0.5851)), 4)
}


//----------------------------------------------------------------------------//
void Tracker::trackit(SP_geometry geo)
{
  Require(geo);
  d_geometry = geo;
  d_X = d_geometry->width_x();
  d_Y = d_geometry->width_y();
  d_Z = d_geometry->width_z();

  // Create database
  d_tracks = new TrackDB(d_quadrature, d_geometry->number_regions());

  // Generate all the start and end points for the tracking procedure
  generate_points();

  // Segmentize all the tracks
  for (size_t a = 0; a < 2 * d_quadrature->number_azimuths_octant(); ++a)
  {
    for (size_t t = 0; t < d_tracks->number_tracks_angle(a); ++t)
    {
      segmentize(d_tracks->track(a, t));
    }
  }
}

//----------------------------------------------------------------------------//
// 2-D Cartesian
void Tracker::generate_points()
{
  using detran_utilities::vec_scale;

  typedef detran_angle::QuadratureFactory   QF;
  typedef QF::SP_basequadrature             SP_squad;

  // spatial quadrature
  SP_squad q;
  d_quadrature->display();
  // do for first quadrant and reflect to get second quadrant points
  size_t na = d_quadrature->number_azimuths_octant();
  for (size_t a = 0; a < na; ++a)
  {
    double phi     = d_quadrature->phi(a);
    double cos_phi = d_quadrature->cos_phi(a);
    double sin_phi = d_quadrature->sin_phi(a);
    double tan_phi = sin_phi / cos_phi;
    double cot_phi = cos_phi / sin_phi;

    size_t dim[] = { 0, 1 };        // the "full" dimension
    double width[] = { d_X, d_Y };  // bounding box

    // By default, the "full" surface is the bottom surface.  If we're too
    // askew, the vertical surface becomes the full surface and we flip
    // all the values.
    double tan_cot_phi = tan_phi;
    double sin_cos_phi = sin_phi;
    double cos_sin_phi = cos_phi;
    if (tan_cot_phi > d_Y / d_X)
    {
      dim[0] = 1;
      dim[1] = 0;
      width[0] = d_Y;
      width[1] = d_X;
      tan_cot_phi = cot_phi;
      sin_cos_phi = cos_phi;
      cos_sin_phi = sin_phi;
    }

    // Primary surface points
    size_t n_P = number_points(phi, width[0], dim[0]);
    q = QF::build_base(d_spatial_quad_type, n_P, 0.0, width[0]);
    vec_dbl P = q->get_x();
    vec_dbl w_P = q->get_w();
    COUT("done primary")

    // Secondary points connected to primary points
    double width_S_p = width[0] * tan_cot_phi;
    q = QF::build_base(d_spatial_quad_type, n_P, 0.0, width_S_p);
    vec_dbl S_p = q->get_x();
    vec_dbl w_S_p = q->get_w();
    COUT("done secondary to primary")

    // Secondary points connected to other secondary points
    double width_S_s = width[1] - width[0] * tan_cot_phi;
    size_t n_S = number_points(phi, width_S_s, dim[1]);
    if (detran_utilities::soft_equiv(width_S_s, 0.0))
    {
      n_S = 0;
    }
    COUT("width_S_s=" << width_S_s << " w1=" << width[1] << " w0=" << width[0])
    vec_dbl S_s, w_S_s;
    if (n_S > 0)
    {
      q = QF::build_base(d_spatial_quad_type, n_S, 0.0, width_S_s);
      S_s = q->get_x();
      w_S_s = q->get_w();
    }
    COUT("done secondary to secondary")

    // Reorder the Y points (and weights) from top to bottom
    if (dim[0] == 0)
    {
      std::reverse(S_p.begin(),   S_p.end());
      std::reverse(w_S_p.begin(), w_S_p.end());
      std::reverse(S_s.begin(),   S_s.end());
      std::reverse(w_S_s.begin(), w_S_s.end());
    }
    else
    {
      std::reverse(P.begin(), P.end());
      std::reverse(w_P.begin(), w_P.end());
    }

    // Add the generated points
    double R0[2], R1[2];
    for (size_t i = 0; i < P.size(); ++i)
    {
      R0[dim[0]] = P[i];                // start: [x, 0] or [0, y]
      R0[dim[1]] = 0.0;
      R1[dim[0]] = width[0];            // end:   [W, y] or [x, H]
      R1[dim[1]] = S_p[i];
      double wt = sin_cos_phi * w_P[i];
      SP_track track(new Track(Point(R0[0], R0[1]), Point(R1[0], R1[1]), wt));
      COUT(*track)
      d_tracks->add_track(a, track);
    }
    COUT("Added primary")
    for (size_t i = 0; i < S_s.size(); ++i)
    {
      R0[dim[0]] = 0.0;
      R0[dim[1]] = S_s[i];              // start: [0, y] or [x, 0]
      R1[dim[0]] = width[0];
      R1[dim[1]] = S_s[i] + width_S_p;  // end:   [W, y+shift] or [x+shift, H]
      double wt = cos_sin_phi * w_S_s[i];
      SP_track track(new Track(Point(R0[0], R0[1]), Point(R1[0], R1[1]), wt));
      d_tracks->add_track(a, track);
    }
    COUT("Added second to second")
    for (size_t i = 0; i < P.size(); ++i)
    {
      R0[dim[0]] = 0.0;
      R0[dim[1]] = S_p[i] + width_S_s;   // start: [0, y] or [x, 0]
      R1[dim[0]] = P[i];
      R1[dim[1]] = width[0];             // end:   [W, y+shift] or [x+shift, H]
      double wt = sin_cos_phi * w_P[i];
      SP_track track(new Track(Point(R0[0], R0[1]), Point(R1[0], R1[1]), wt));
      d_tracks->add_track(a, track);
    }
    COUT("Added second to primary")
  }

  // Mirror the tracks in the second quadrant.
  for (size_t a = 0; a < na; ++a)
  {
    for (size_t t = 0; t < d_tracks->number_tracks_angle(a); ++t)
    {
      SP_track t0 = d_tracks->track(a, t);
      Point P0 = Point(d_X - t0->enter().x(), t0->enter().y());
      Point P1 = Point(d_X - t0->exit().x(), t0->exit().y());
      SP_track t1(new Track(P0, P1, t0->width()));
      d_tracks->add_track(na + a, t1);

    }
  }
  //d_tracks->display();
}

//----------------------------------------------------------------------------//
// CSG SEGMENTATION SUPPORT
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// contain a segment with its midpoint along the track
struct MidpointSegment
{
  MidpointSegment(const Point &p, const size_t reg, const double length)
    : midpoint(p), segment(reg, length) {/* ... */}
  Point   midpoint;
  Segment segment;
};

//----------------------------------------------------------------------------//
// compares two segments based on their midpoints w/r to the track origin
struct MidpointSegmentCompare
{
  MidpointSegmentCompare(const Point &o) : d_origin(o) {/* ... */}
  bool operator() (const MidpointSegment &r0, const MidpointSegment &r1)
  {
    return distance(r0.midpoint, d_origin) < distance(r1.midpoint, d_origin);
  }
  const Point &d_origin;
};

//----------------------------------------------------------------------------//
void Tracker::segmentize(SP_track track)
{

  // starting point, end point, track length, and direction
  Point P0 = track->enter();
  Point P1 = track->exit();
  double L = distance(P0, P1);
  Point D  = (P1 - P0) * (1.0 / L);

  // vector of segments with their midpoints
  std::vector<MidpointSegment> segments;

  // loop through all the regions
  for (size_t r = 0; r < d_geometry->number_regions(); ++r)
  {
    //COUT("region: " << r)
    Region::SP_region region = d_geometry->region(r);

    // get all intersections with the node constituents.  note, for simplicity,
    // we start the ray behind the starting point to catch it as an intercept
    CSG_Node::vec_point points =
      region->top_node()->intersections(P0 - 1.0e-8*D, D, L);

    if (points.size() == 0) continue;

    // loop through all points (skipping the pathological case of one point)
    for (size_t i = 0; i < points.size() - 1; ++i)
    {

//      COUT(points[i] << points[i+1] <<
//           " " << distance(points[i], points[i+1]) <<
//           " " << distance(points[i  ], P0) <<
//           " " << distance(points[i+1], P0))
//
      if (abs(distance(points[i], points[i+1]) < 1e-14)) continue;

      // make sure we have proper order
      Assert(distance(points[i], P0) < distance(points[i+1], P0));

      // compute the midpoint
      Point midpoint = 0.5 * (points[i] + points[i+1]);

      // compute and add the segment if the region contains the midpoint
      if (region->top_node()->contains(midpoint))
      {
        double length = distance(points[i], points[i+1]);
        //COUT("length = " << length)
        segments.push_back(MidpointSegment(midpoint, r, length));
      }

    } // end points

  } // end regions

  // sort the segments in the order they appear in the track
  std::sort(segments.begin(), segments.end(), MidpointSegmentCompare(P0));

  using detran_utilities::soft_equiv;
  double track_length = 0;
  double ref_track_length = distance(track->enter(), track->exit());
  for (size_t i = 0; i < segments.size(); ++i)
    track_length += segments[i].segment.length();

  Ensurev(soft_equiv(track_length, ref_track_length, 1e-12),
          "Actual length:   " + AsString(ref_track_length) +
          "Computed length: " + AsString(track_length));

  // add them to the track
  for (size_t i = 0; i < segments.size(); ++i)
    track->add_segment(segments[i].segment);
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Tracker.cc
//----------------------------------------------------------------------------//
