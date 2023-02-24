//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Tracker.cc
 *  @brief Tracker class member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Tracker.hh"
#include "QuadraticSurfaceFactory.hh"
#include "angle/QuadratureFactory.hh"
#include "utilities/SoftEquivalence.hh"
#include "utilities/MathUtilities.hh"
#include <iostream>
#include <cmath>
#include <cfloat>

#define TRKRMSG(c) std::cout << c << std::endl;

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Tracker::Tracker(SP_db db, SP_quadrature q)
  : d_quadrature(q)
  , d_X(0)
  , d_Y(0)
  , d_Z(0)
  , d_average_width(0.1)
  , d_minimum_number(1)
  , d_spatial_quad_type("uniform")
{
  Require(db);
  Require(d_quadrature);
  if (db->check("tracker_maximum_spacing"))
  {
    d_average_width = db->get<double>("tracker_maximum_spacing");
    Assert(d_average_width > 0.0);
  }
  if (db->check("tracker_minimum_number"))
  {
    d_minimum_number = db->get<int>("tracker_minimum_number");
    Assert(d_minimum_number > 0);
  }
  if (db->check("tracker_spatial_quad_type"))
  {
    d_spatial_quad_type = db->get<std::string>("tracker_spatial_quad_type");
  }
  if (db->check("tracker_symmetric_tracks"))
  {
    d_symmetric_tracks = 0 != db->get<int>("tracker_symmetric_tracks");
  }
  d_trackdb = std::make_shared<TrackDB>(q);
}

//----------------------------------------------------------------------------//
// IMPLEMENTATION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
std::shared_ptr<TrackDB> Tracker::trackit(std::shared_ptr<Geometry> geo)
{
  Require(geo);
  d_geometry = geo;
  d_X = d_geometry->width_x();
  d_Y = d_geometry->width_y();
  d_Z = d_geometry->width_z();

  TRKRMSG("generating track points...")
  if (d_quadrature->dimension() == 2)
    generate_tracks_2D_cartesian();
  else
    generate_tracks_3D_cartesian();
  TRKRMSG("segmentizing...")
  for (auto& angle_track: d_trackdb->tracks())
    for (auto &track: angle_track.second)
      segmentize(track);
  TRKRMSG("...tracking complete.");

  std::shared_ptr<TrackDB> track_db;
  return track_db;
}

//----------------------------------------------------------------------------//
std::shared_ptr<TrackDB> Tracker::trackit(Mesh::SP_mesh mesh)
{
  auto geo = std::make_shared<Geometry>(mesh);
  return trackit(geo);
}

//----------------------------------------------------------------------------//
Point Tracker::find_exit(const Point &entry, const Point &direction)
{
  Require(direction.x() > 0.0);
  Require(direction.y() > 0.0);
  Require(direction.z() >= 0.0);

  auto X_plane_east = QuadraticSurfaceFactory::CreatePlaneX(d_X);
  auto Y_plane_north = QuadraticSurfaceFactory::CreatePlaneY(d_Y);
  auto Z_plane_top = QuadraticSurfaceFactory::CreatePlaneY(d_Z);

  Ray ray(entry, direction);

  auto intersection_X = X_plane_east->intersections(ray);
  Ensure(intersection_X.size() == 1);
  Point exit = intersection_X[0];
  double dist = distance(intersection_X[0], entry);

  auto intersection_Y = Y_plane_north->intersections(ray);
  Ensure(intersection_Y.size() == 1);
  double dist_Y = distance(intersection_Y[0], entry);
  if (dist_Y < dist)
  {
    dist = dist_Y;
    exit = intersection_Y[0];
  }
  if (direction.z() != 0)
  {
    auto intersection_Z = Z_plane_top->intersections(ray);
    Ensure(intersection_Z.size() == 1);
    double dist_Z = distance(intersection_Z[0], entry);
    if (dist_Z < dist) exit = intersection_Z[0];
  }
  return exit;
}

//----------------------------------------------------------------------------//
void Tracker::generate_tracks_2D_cartesian()
{
	// For each azimuth in the first quadrant, points are distributed along
	// each incident face according to an n-point quadrature. Here, those
	// n points come from a target average track spacing.  So, for anything but
	// uniform quadrature points, some spacings will be larger.
	//
	// Once the initial points are selected, symmetry is used to find the
	// tracks for the remaining three coordinates.
	typedef detran_angle::QuadratureFactory   QF;
	auto& Q = *d_quadrature;

	auto& d_tracks = d_trackdb->tracks();

	// generate first-quadrant tracks
	for (int a = 0; a < Q.number_azimuths_octant(); ++a)
	{
	  // sanity checks to ensure product quadrature implementation
	  // doesn't break things
	  std::vector<Track::SP_track> tracks;

		// number of points on each face
		int num_horz = std::max(d_minimum_number, (int)std::ceil(d_X / (d_average_width / Q.cos_phi(a))));
		int num_vert = std::max(d_minimum_number, (int)std::ceil(d_Y / (d_average_width / Q.sin_phi(a))));

		// define quadratures
		auto q_horz = QF::build_base(d_spatial_quad_type, num_horz, 0.0, d_X);
		auto q_vert = QF::build_base(d_spatial_quad_type, num_vert, 0.0, d_Y);
		{
		  auto &y = q_vert->get_x();
		  std::vector<double> x(y.size(), 0.0);
		  auto &w = q_vert->get_w();
		  for (int j = 0; j < y.size(); ++j)
		  {
		    Point entry(x[j], y[j]);
		    Point exit = find_exit(entry, Point(Q.cos_phi(a), Q.sin_phi(a)));
		    auto track = std::make_shared<Track>(entry, exit, w[j]*Q.sin_phi(a));
		    tracks.push_back(track);
		  }
		}
		{
		  auto &x = q_horz->get_x();
		  std::vector<double> y(x.size(), 0.0);
		  auto &w = q_horz->get_w();
		  for (int i = 0; i < x.size(); ++i)
		  {
        Point entry(x[i], y[i]);
        Point exit = find_exit(entry, Point(Q.cos_phi(a), Q.sin_phi(a)));
        auto track = std::make_shared<Track>(entry, exit, w[i]*Q.cos_phi(a));
        tracks.push_back(track);
		  }
		}
		d_tracks[{Q.angle(a, 0), 0}] = tracks;
	}

	// create mirrors such that, e.g., entry.x() becomes d_X - entry.x()
	std::map<std::pair<int,int>,std::vector<Track::SP_track>> mirror_tracks;

  for (auto angle_track: d_tracks)
  {
    std::vector<Track::SP_track>& tracks = angle_track.second;
    int a = angle_track.first.first;
    std::vector<Track::SP_track> new_tracks;
    for (auto track: tracks)
    {
      const auto &enter = track->enter();
      auto new_enter = Point(d_X-enter.x(), enter.y());
      const auto &exit = track->exit();
      auto new_exit = Point(d_X-exit.x(), exit.y());
      auto new_track = std::make_shared<Track>(new_enter, new_exit, track->spatial_weight());
      new_tracks.push_back(track);
    }
    Assert(Q.mu(0, a)==-Q.mu(1, a));
    Assert(Q.eta(0, a)==Q.eta(1, a));
    d_tracks[{a, 1}] = new_tracks;
  }

  TRKRMSG(" all tracks initiated for 2-d problem. ")
}


//----------------------------------------------------------------------------//
/**
 *   For 3-D, we use ray tracing on the bounding box.  For each incident
 *   face, a product spatial quadrature is used.  This generally leads
 *   to asymmetric tracking when symmetry is applied to reduce the tracking
 *   info to the first four octants (since a track from 6 to 0 can be used
 *   from 0 to 6).  If symmetric tracking is desired, separate tracks should
 *   be kept for up to down tracks.
 */
void Tracker::generate_tracks_3D_cartesian()
{
//  using detran_utilities::soft_equiv;
//  using detran_utilities::vec_scale;
//  using std::max;
//  using std::min;
//  using std::ceil;
//
//  typedef detran_angle::QuadratureFactory   QF;
//  typedef QF::SP_basequadrature             SP_squad;
//
//  // spatial quadrature
//  SP_squad q;
//
//  // point representing true outer bounds and the max track length to get there
//  Point bound(d_X, d_Y, d_Z);
//  double max_length = distance(bound, Point(0, 0, 0)) + 0.001;
//
//  // loop over all first octant angles
//  for (size_t a = 0; a < d_quadrature->number_azimuths_octant(); ++a)
//  {
//    COUT("a = " << a)
//    for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
//    {
//      COUT("p = " << p)
//
//      // generate points along x, y, and z (assuming even spacing for defining
//      // the number of points to use; for g-l, use larger spacing)
//      size_t angle = d_quadrature->angle(a, p);
//      size_t num_points[3]; // number of entrance points along each axis
//      vec_dbl p_xyz[3];     // the x, y, and z points
//      vec_dbl w_xyz[3];     //   and corresponding weights
//      for (size_t d = 0; d < 3; ++d)
//      {
//        double sin_xyz = sqrt(1.0 - pow(d_quadrature->cosines(d)[angle], 2));
//        num_points[d] = max((int)ceil(bound[d]/d_maximum_spacing/sin_xyz), 3);
//        q = QF::build_base(d_spatial_quad_type, num_points[d], 0.0, bound[d]);
//        p_xyz[d] = q->get_x();
//        w_xyz[d] = q->get_w();
//      }
//
//      // track direction
//      Point D(d_quadrature->mu(0, angle),
//              d_quadrature->eta(0, angle),
//              d_quadrature->xi(0, angle));
//
//      // loop over the two dimensions of the incident surface
//      size_t remaining_d[3] = {2, 1, 0};
//      for (size_t d0 = 0; d0 < 3; ++d0)
//      {
//        for (size_t d1 = d0+1; d1 < 3; ++d1)
//        {
//          // determine the dimension perpendicular to this surface
//          size_t d2 = remaining_d[d0+d1-1];
//          for (size_t i = 0; i < p_xyz[d0].size(); ++i)
//          {
//            for (size_t j = 0; j < p_xyz[d1].size(); ++j)
//            {
//              // entrance point
//              Point P0;
//              P0[d0] = p_xyz[d0][i];
//              P0[d1] = p_xyz[d1][j];
//              // exit point
//              double t = DBL_MAX;
//              for (size_t s = 0; s < 3; ++s)
//                t = std::min(t, (bound[s] - P0[s]) / D[s]);
//              Assert(!soft_equiv(t, DBL_MAX) && t > 0.0);
//              Point P1 = P0 + t * D;
//              // track spatial weight = cross sectional "area"
//              double wt = w_xyz[d0][i] * w_xyz[d1][j];
//              d_tracks->add_track(a, p, SP_track(new Track(P0, P1, wt)));
//            }
//          }
//        }
//      }
//    } // end polar
//  } // end azimuth
//
//  // apply symmetry (note, this must be redone if backward angles are to
//  // have the same entrance point distribution as the forward angles)
//  size_t na = d_quadrature->number_azimuths_octant();
//  for (size_t a = 0; a < na; ++a)
//  {
//    for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
//    {
//      for (size_t t = 0; t < d_tracks->number_tracks(a); ++t)
//      {
//        SP_track t0 = d_tracks->track(a, p, t);
//        Point P0 = t0->enter(), P1 = t0->exit();
//        for (size_t q = 1; q < 4; ++q)
//        {
//          if (q == 1 || q == 2)
//          {
//            P0[0] = d_X - P0[0];
//            P1[0] = d_X - P1[0];
//          }
//          if (q == 2 || q == 3)
//          {
//            P0[1] = d_Y - P0[1];
//            P1[1] = d_Y - P1[1];
//          }
//          SP_track t1(new Track(P0, P1, t0->width()));
//          d_tracks->add_track(a + q * na, p, t1);
//        }
//      }
//    }
//  }

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
void Tracker::segmentize(std::shared_ptr<Track> track)
{
  using detran_utilities::soft_equiv;
  typedef CSG_Node::vec_point vec_point;

  // starting point, end point, track length, and direction
  Point P0 = track->enter();
  Point P1 = track->exit();
  TRKRMSG("enter=" << P0 << " exit=" << P1)
  double L = distance(P0, P1);
  Point  D = (P1 - P0) * (1.0 / L);

  // define the ray to start just slightly behind the entrance so the first
  // intersection is automatically included in the intersections, and correct
  // the maximum length to go just beyond the exit to ensure it's included.
  double eps = 1.0e-9;
  Ray ray(P0 - eps*D, D);
  double ray_length = L + 2.0*eps;

  // vector of segments with their midpoints
  std::vector<MidpointSegment> segments;

  // loop through all the regions
  bool found = false;
  for (size_t r = 0; r < d_geometry->number_regions(); ++r)
  {
    auto region = d_geometry->region(r);

    if (!region->intersects_bounding_box(ray, ray_length))
    {
      TRKRMSG("I'm a ducking piece of ship.  Puck this box:"
            << region->bound_min() << region->bound_max())
      continue;
    }
    TRKRMSG("In BBox so need points"
          << " P0="   << P0
          << " D=" << D
          << " ray = " << ray.origin << " " << ray.direction
          << " min=" << region->bound_min()
          << " max=" << region->bound_max())
    TRKRMSG(" ")
    // get all intersections with the node constituents
    vec_point points = region->top_node()->intersections(ray, ray_length);

    Assert(points.size() > 0);
    found = true;
    // loop through all points (skipping the pathological case of one point)
    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      TRKRMSG(points[i])
      if (soft_equiv(std::abs(distance(points[i], points[i+1])), 0.0)) continue;

      // make sure we have proper order
      Assert(distance(points[i], P0) < distance(points[i+1], P0));

      // compute the midpoint
      Point midpoint = 0.5 * (points[i] + points[i+1]);

      // compute and add the segment if the region contains the midpoint
      if (region->contains(midpoint))
      {
        double segment_length = distance(points[i], points[i+1]);
        segments.push_back(MidpointSegment(midpoint, r, segment_length));
      }
    } // end points
  } // end regions
  Assert(found);
  // sort the segments in the order they appear in the track
  std::sort(segments.begin(), segments.end(), MidpointSegmentCompare(P0));

#ifdef DETRAN_ENABLE_DEBUG
  using detran_utilities::soft_equiv;
  double track_length = 0;
  double ref_track_length = distance(track->enter(), track->exit());
  for (size_t i = 0; i < segments.size(); ++i)
    track_length += segments[i].segment.length();
  Ensurev(soft_equiv(track_length, ref_track_length, 1e-12),
          "Actual length:   " + AsString(ref_track_length) +
          "Computed length: " + AsString(track_length));
#endif

  // add them to the track
  for (size_t i = 0; i < segments.size(); ++i)
    track->add_segment(segments[i].segment);
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Tracker.cc
//----------------------------------------------------------------------------//
