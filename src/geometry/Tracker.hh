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
#include <memory>
#include <vector>
#include <map>

namespace detran_geometry
{

/**
 *  @class Tracker
 *  @brief Track a Cartesian mesh or a CSG-based geometry in two
 *         or three dimensions.
 *
 *  The Tracker is responsible for generating tracks across the given
 *  geometry.
 *
 *  Relevant database enteries:
 *    - tracker_maximum_spacing   [double]
 *    - tracker_spatial_quad_type [string]
 *    - tracker_normalize_lengths [int]
 *    - tracker_symmetric_tracks  [int]
 */
class GEOMETRY_EXPORT Tracker
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input             SP_db;
  typedef detran_angle::ProductQuadrature::SP_quadrature  SP_quadrature;
  typedef std::vector<std::shared_ptr<Track>>                 vec_track;
  typedef std::map<std::pair<int, int>, vec_track>            map_track;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    db    Parameter database
   *  @param    q     Product quadrature used for angular distribution of tracks
   */
  Tracker(SP_db db, SP_quadrature q);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Track a Cartesian mesh
  std::shared_ptr<TrackDB> trackit(std::shared_ptr<Mesh> mesh);

  /// Track a CSG tree with a Cartesian bounding box
  std::shared_ptr<TrackDB> trackit(std::shared_ptr<Geometry> geo);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Angular quadrature that *must* be a product quadrature
  SP_quadrature d_quadrature;
  /// Desired average distance between any two parallel tracks
  double d_average_width;
  /// Minimum number of entry points on an incident surface for any angle
  int d_minimum_number;
  /// Quadrature used to define track spatial distribution
  std::string d_spatial_quad_type;
  /// Flag to indicate symmetric tracking is to be done in 3-D
  bool d_symmetric_tracks;
  /// Geometry to be tracked
  std::shared_ptr<Geometry> d_geometry;
  //@{
  ///  Bounding box dimensions
  double d_X;
  double d_Y;
  double d_Z;
  //@}

  /// map of octant,angle indices to vector of tracks.
  std::shared_ptr<TrackDB> d_trackdb;
  map_track d_tracks;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  // Given an entry and direction on the bounding box, find the exit.
  Point find_exit(const Point &entry, const Point &direction);

  /// generate track points for a 2D cartesian geometry
  void generate_tracks_2D_cartesian();

  /// generate track points for a 3D cartesian geometry
  void generate_tracks_3D_cartesian();

  /// cast a track across the domain and segmentize it
  void segmentize(Track::SP_track track);

};

} // end namespace detran_geometry

#endif /* detran_geometry_TRACKER_HH_ */

//----------------------------------------------------------------------------//
//              end of file Tracker.hh
//----------------------------------------------------------------------------//
