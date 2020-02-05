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
  if (d_db->check("tracker_symmetric_tracks"))
  {
    d_symmetric_tracks = 0 != d_db->get<int>("tracker_symmetric_tracks");
  }
}

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

  size_t n = std::max((int)tmp, 3); // minimum of 3 tracks per surface
  return n;

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
  d_tracks = new TrackDB(d_quadrature);


  if (d_quadrature->dimension() == 2)
  {
    // Generate all track entrances and exits
    generate_points_2D_cartesian();

    // Sort the points.
    d_tracks->sort();

    // Segmentize all the tracks
    for (size_t a = 0; a < 2 * d_quadrature->number_azimuths_octant(); ++a)
      for (size_t t = 0; t < d_tracks->number_tracks(a); ++t)
        segmentize(d_tracks->track(a, 0, t));
  }
  else
  {
    // Generate all track entrances and exits
    generate_points_3D_cartesian();

    // Sort the points.
    d_tracks->sort();

    // Segmentize all the tracks
    for (size_t a = 0; a < 4 * d_quadrature->number_azimuths_octant(); ++a)
      for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
        for (size_t t = 0; t < d_tracks->number_tracks(a); ++t)
          segmentize(d_tracks->track(a, 0, t));
  }
}

//----------------------------------------------------------------------------//
void Tracker::trackit(SP_mesh mesh)
{
  SP_geometry geo = mesh_to_geometry(mesh);
  trackit(geo);
}

Tracker::SP_geometry Tracker::mesh_to_geometry(SP_mesh mesh)
{
  typedef QuadraticSurfaceFactory QSF;

  Require(mesh);
  Require(mesh->dimension() > 1);

  d_mesh = mesh;
  SP_geometry geo(new Geometry(d_mesh->total_width_x(),
                               d_mesh->total_width_y(),
                               d_mesh->total_width_z()));

  // create a plane-based geometry out of the mesh
  vec_dbl edges[3];
  Surface::vec_surface surfaces[3];
  for (size_t d = 0; d < 3; ++d)
  {
    size_t n = d_mesh->number_cells(d);
    edges[d].resize(n + 1);
    surfaces[d].resize(n + 1);
    size_t c[3] = {0, 0, 0};
    c[d] = 1;
    double edge = 0.0;
    for (size_t i = 0; i <= n; ++i)
    {
      edges[d][i]    = edge;
      surfaces[d][i] = QSF::CreatePlane(c[0], c[1], c[2], edge);
      if (i < n) edge += d_mesh->width(d, i);
    }
  }

  const vec_int &mat_map = d_mesh->mesh_map("MATERIAL");
  for (size_t k = 0; k < d_mesh->number_cells_z(); ++k)
  {
    for (size_t j = 0; j < d_mesh->number_cells_y(); ++j)
    {
      for (size_t i = 0; i < d_mesh->number_cells_x(); ++i)
      {
        size_t m    = mat_map[d_mesh->index(i, j, k)];
        Point b_min, b_max;
        b_min = Point(edges[0][i  ], edges[1][j  ], edges[2][k  ]) - 1e-5;
        b_max = Point(edges[0][i+1], edges[1][j+1], edges[2][k+1]) + 1e-5;
        SP_region r = Region::Create(m, b_min, b_max);
        r->append(surfaces[0][i  ], true);
        r->append(surfaces[0][i+1], false);
        r->append(surfaces[1][j  ], true);
        r->append(surfaces[1][j+1], false);
        if (d_mesh->dimension() == 3) r->append(surfaces[2][k  ], true);
        if (d_mesh->dimension() == 3) r->append(surfaces[2][k+1], false);
        geo->add_region(r);
      }
    }
  }

  return geo;
}


//----------------------------------------------------------------------------//
void Tracker::generate_points_2D_cartesian()
{
  using detran_utilities::vec_scale;

  typedef detran_angle::QuadratureFactory   QF;
  typedef QF::SP_basequadrature             SP_squad;

  // spatial quadrature
  SP_squad q;

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

    // Secondary points connected to primary points
    double width_S_p = width[0] * tan_cot_phi;
    q = QF::build_base(d_spatial_quad_type, n_P, 0.0, width_S_p);
    vec_dbl S_p = q->get_x();
    vec_dbl w_S_p = q->get_w();

    // Secondary points connected to other secondary points
    double width_S_s = width[1] - width[0] * tan_cot_phi;
    size_t n_S = number_points(phi, width_S_s, dim[1]);
    if (detran_utilities::soft_equiv(width_S_s, 0.0))
    {
      n_S = 0;
    }

    vec_dbl S_s, w_S_s;
    if (n_S > 0)
    {
      q = QF::build_base(d_spatial_quad_type, n_S, 0.0, width_S_s);
      S_s = q->get_x();
      w_S_s = q->get_w();
    }

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
      d_tracks->add_track(a, 0, track);
    }
    for (size_t i = 0; i < S_s.size(); ++i)
    {
      R0[dim[0]] = 0.0;
      R0[dim[1]] = S_s[i];              // start: [0, y] or [x, 0]
      R1[dim[0]] = width[0];
      R1[dim[1]] = S_s[i] + width_S_p;  // end:   [W, y+shift] or [x+shift, H]
      double wt = cos_sin_phi * w_S_s[i];
      SP_track track(new Track(Point(R0[0], R0[1]), Point(R1[0], R1[1]), wt));
      d_tracks->add_track(a, 0, track);
    }
    for (size_t i = 0; i < P.size(); ++i)
    {
      R0[dim[0]] = 0.0;
      R0[dim[1]] = S_p[i] + width_S_s;   // start: [0, y] or [x, 0]
      R1[dim[0]] = P[i];
      R1[dim[1]] = width[0];             // end:   [W, y+shift] or [x+shift, H]
      double wt = sin_cos_phi * w_P[i];
      SP_track track(new Track(Point(R0[0], R0[1]), Point(R1[0], R1[1]), wt));
      d_tracks->add_track(a, 0, track);
    }
  }

  // Mirror the tracks in the second quadrant.
  for (size_t a = 0; a < na; ++a)
  {
    for (size_t t = 0; t < d_tracks->number_tracks(a); ++t)
    {
      SP_track t0 = d_tracks->track(a, 0, t);
      Point P0 = Point(d_X - t0->enter().x(), t0->enter().y());
      Point P1 = Point(d_X - t0->exit().x(), t0->exit().y());
      SP_track t1(new Track(P0, P1, t0->width()));
      d_tracks->add_track(na + a, 0, t1);

    }
  }

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
void Tracker::generate_points_3D_cartesian()
{
  using detran_utilities::soft_equiv;
  using detran_utilities::vec_scale;
  using std::max;
  using std::min;
  using std::ceil;

  typedef detran_angle::QuadratureFactory   QF;
  typedef QF::SP_basequadrature             SP_squad;

  // spatial quadrature
  SP_squad q;

  // point representing true outer bounds and the max track length to get there
  Point bound(d_X, d_Y, d_Z);
  double max_length = distance(bound, Point(0, 0, 0)) + 0.001;

  // loop over all first octant angles
  for (size_t a = 0; a < d_quadrature->number_azimuths_octant(); ++a)
  {
    COUT("a = " << a)
    for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
    {
      COUT("p = " << p)

      // generate points along x, y, and z (assuming even spacing for defining
      // the number of points to use; for g-l, use larger spacing)
      size_t angle = d_quadrature->angle(a, p);
      size_t num_points[3]; // number of entrance points along each axis
      vec_dbl p_xyz[3];     // the x, y, and z points
      vec_dbl w_xyz[3];     //   and corresponding weights
      for (size_t d = 0; d < 3; ++d)
      {
        double sin_xyz = sqrt(1.0 - pow(d_quadrature->cosines(d)[angle], 2));
        num_points[d] = max((int)ceil(bound[d]/d_maximum_spacing/sin_xyz), 3);
        q = QF::build_base(d_spatial_quad_type, num_points[d], 0.0, bound[d]);
        p_xyz[d] = q->get_x();
        w_xyz[d] = q->get_w();
      }

      // track direction
      Point D(d_quadrature->mu(0, angle),
              d_quadrature->eta(0, angle),
              d_quadrature->xi(0, angle));

      // loop over the two dimensions of the incident surface
      size_t remaining_d[3] = {2, 1, 0};
      for (size_t d0 = 0; d0 < 3; ++d0)
      {
        for (size_t d1 = d0+1; d1 < 3; ++d1)
        {
          // determine the dimension perpendicular to this surface
          size_t d2 = remaining_d[d0+d1-1];
          for (size_t i = 0; i < p_xyz[d0].size(); ++i)
          {
            for (size_t j = 0; j < p_xyz[d1].size(); ++j)
            {
              // entrance point
              Point P0;
              P0[d0] = p_xyz[d0][i];
              P0[d1] = p_xyz[d1][j];
              // exit point
              double t = DBL_MAX;
              for (size_t s = 0; s < 3; ++s)
                t = std::min(t, (bound[s] - P0[s]) / D[s]);
              Assert(!soft_equiv(t, DBL_MAX) && t > 0.0);
              Point P1 = P0 + t * D;
              // track spatial weight = cross sectional "area"
              double wt = w_xyz[d0][i] * w_xyz[d1][j];
              d_tracks->add_track(a, p, SP_track(new Track(P0, P1, wt)));
            }
          }
        }
      }
    } // end polar
  } // end azimuth

  // apply symmetry (note, this must be redone if backward angles are to
  // have the same entrance point distribution as the forward angles)
  size_t na = d_quadrature->number_azimuths_octant();
  for (size_t a = 0; a < na; ++a)
  {
    for (size_t p = 0; p < d_quadrature->number_polar_octant(); ++p)
    {
      for (size_t t = 0; t < d_tracks->number_tracks(a); ++t)
      {
        SP_track t0 = d_tracks->track(a, p, t);
        Point P0 = t0->enter(), P1 = t0->exit();
        for (size_t q = 1; q < 4; ++q)
        {
          if (q == 1 || q == 2)
          {
            P0[0] = d_X - P0[0];
            P1[0] = d_X - P1[0];
          }
          if (q == 2 || q == 3)
          {
            P0[1] = d_Y - P0[1];
            P1[1] = d_Y - P1[1];
          }
          SP_track t1(new Track(P0, P1, t0->width()));
          d_tracks->add_track(a + q * na, p, t1);
        }
      }
    }
  }

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
  using detran_utilities::soft_equiv;
  typedef CSG_Node::vec_point vec_point;

  // starting point, end point, track length, and direction
  Point P0 = track->enter();
  Point P1 = track->exit();
  COUT("enter=" << P0 << " exit=" << P1)
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
    Region::SP_region region = d_geometry->region(r);

    if (!region->intersects_bounding_box(ray, ray_length))
    {
      COUT("I'm a fucking piece of shit.  Fuck this box:"
            << region->bound_min() << region->bound_max())
      continue;
    }
    COUT("In BBox so need points"
          << " P0="   << P0
          << " D=" << D
          << " ray = " << ray.origin << " " << ray.direction
          << " min=" << region->bound_min()
          << " max=" << region->bound_max())
    COUT(" ")
    // get all intersections with the node constituents
    vec_point points = region->top_node()->intersections(ray, ray_length);

    Assert(points.size() > 0);
    found = true;
    // loop through all points (skipping the pathological case of one point)
    for (size_t i = 0; i < points.size() - 1; ++i)
    {
      COUT(points[i])
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
