//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Tracker.hh
 *  @brief Tracker class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_TRACKER_HH_
#define detran_geometry_TRACKER_HH_

#include "Track.hh"
#include "TrackDB.hh"
#include "Mesh.hh"
#include "Geometry.hh"
#include "angle/ProductQuadrature.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_geometry
{

/**
 *  @class Tracker
 *  @brief Track a Cartesian mesh or a CSG-based geometry
 *
 *  Relevant database enteries:
 *    - tracker_maximum_spacing   [double]
 *    - tracker_spatial_quad_type [string]
 *    - tracker_normalize_lengths [int]
 *    - tracker_symmetric_tracks  [int]
 *
 */
class GEOMETRY_EXPORT Tracker
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Tracker>                   SP_tracker;
  typedef detran_utilities::InputDB::SP_input             SP_db;
  typedef Mesh::SP_mesh                                   SP_mesh;
  typedef Geometry::SP_geometry                           SP_geometry;
  typedef detran_angle::ProductQuadrature::SP_quadrature  SP_quadrature;
  typedef Track::SP_track                                 SP_track;
  typedef TrackDB::SP_trackdb                             SP_trackdb;
  typedef detran_utilities::vec_dbl                       vec_dbl;
  typedef detran_utilities::vec2_dbl                      vec2_dbl;
  typedef detran_utilities::vec_int                       vec_int;
  typedef detran_utilities::size_t                        size_t;
  typedef Region::SP_region                               SP_region;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    db    Parameter database
   *  @param    q     Product quadrature
   */
  Tracker(SP_db db, SP_quadrature q);

  /// SP constructor
  static SP_tracker Create(SP_db db, SP_quadrature q);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Track a Cartesian mesh
  void trackit(SP_mesh mesh);

  /// Track a CSG tree
  void trackit(SP_geometry geo);

  SP_trackdb trackdb() const
  {
    Require(d_tracks);
    return d_tracks;
  }

  // Normalize the track segments based on actual volumes.
  void normalize();


  double maximum_spacing() const {return d_maximum_spacing;}

  SP_geometry mesh_to_geometry(SP_mesh mesh);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Parameter database
  SP_db d_db;
  /// Angular quadrature that *must* be a product quadrature
  SP_quadrature d_quadrature;
  /// Maximum desired distance between any two parallel tracks
  double d_maximum_spacing;
  /// Quadrature used to define track spatial distribution
  std::string d_spatial_quad_type;
  /// Flag to indicate symmetric tracking is to be done in 3-D
  bool d_symmetric_tracks;
  /// Mesh to be tracked
  SP_mesh d_mesh;
  /// Geometry to be tracked
  SP_geometry d_geometry;
  //@{
  ///  Bounding box dimensions
  double d_X;
  double d_Y;
  double d_Z;
  //@}
  /// Database of tracks
  SP_trackdb d_tracks;
  // Number azimuths per octant
  vec_dbl d_x;
  vec_dbl d_y;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  // compute number of points on a (segment of one) side
  size_t number_points(const double phi,
                       const double L,
                       const size_t flag = 0);

  /// generate track points for a 2D cartesian geometry
  void generate_points_2D_cartesian();

  /// generate track points for a 3D cartesian geometry
  void generate_points_3D_cartesian();

  /// cast a track across the domain and segmentize it
  void segmentize(SP_track track);

};

} // end namespace detran_geometry

#endif /* detran_geometry_TRACKER_HH_ */

//----------------------------------------------------------------------------//
//              end of file Tracker.hh
//----------------------------------------------------------------------------//
