//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TrackDB.hh
 *  @brief TrackDB class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_TRACKDB_HH_
#define detran_geometry_TRACKDB_HH_

#include "geometry/geometry_export.hh"
#include "geometry/Track.hh"
#include "angle/ProductQuadrature.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_geometry
{

/**
 *  @class TrackDB
 *  @brief Database of tracks.
 *
 *  The track database contains all the information needed to represent the
 *  MOC problem with respect to space and angle.
 *
 *  Currently, this is limited to 2-D tracking.  Tracks are only assigned
 *  for azimuths over [0, pi], using symmetry for [pi, 2*pi].
 *
 *  Tracks are ordered in space from left to right w/r to the incident
 *  surface
 *
 *  Tracks are stored in order in order of decreasing
 *  y and increasing x for 0, pi/2 and decreasing
 *  x for pi/2, pi.  This is how tracks are indexed.
 *  The track index in the other two octancts keeps
 *  the index of their reflection.
 *
 */
/**
 *  @example geometry/test/test_TrackDB.cc
 *
 *  Test of TrackDB class
 */
class GEOMETRY_EXPORT TrackDB
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<TrackDB>                   SP_trackdb;
  typedef detran_angle::ProductQuadrature::SP_quadrature  SP_quadrature;
  typedef Track::SP_track                                 SP_track;
  typedef std::vector<SP_track>                           vec_track;
  typedef std::vector<vec_track>                          vec2_track;
  typedef detran_utilities::size_t                        size_t;
  typedef detran_utilities::vec_int                       vec_int;
  typedef detran_utilities::vec_dbl                       vec_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    quad           product quadrature
   *  @param    num_regions    number of flat source regions
   */
  TrackDB(SP_quadrature quad,
          const size_t  num_regions);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Get a track
   *  @param    a       azimuth index
   *  @param    t       track index within azimuths
   */
  SP_track track(const size_t a, const size_t t);

  /**
   *  @brief Number of tracks in an azimuthal angle
   *  @param    a       azimuth index
   */
  size_t number_tracks_angle(const size_t a) const;

  /**
   *  @brief Total number azimuths tracked
   *  @param    a       azimuth index
   *  @param    t       track index within azimuths
   */
  size_t number_angles() const;

  /**
   *  @brief Add a track
   *  @param    a       azimuth index
   *  @param    t       track index within azimuths
   */
  void add_track(const size_t a, SP_track t);

  /**
   *  @brief Define information for an azimuth
   *  @param    a       azimuth index
   *  @param    t       track index within azimuths
   */
  void setup_angle(const size_t a,
                   const double c_phi,
                   const double s_phi);

  /// Normalize the tracks given a vector of true volumes.
  void normalize(const vec_dbl &volume);

  /// Pretty display of all track
  void display() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Quadrature
  SP_quadrature d_quadrature;
  /// Number of azimuths in first two octants
  size_t d_number_azimuths;
  /// Number of flat source regions
  size_t d_number_regions;
  /// Tracks by [azimuth][space]
  vec2_track d_tracks;
  /// Azimuthal cosines.
  vec_dbl d_cos_phi;
  /// Azimuthal sines.
  vec_dbl d_sin_phi;

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<TrackDB>)

} // end namespace detran_geometry

#endif // detran_geometry_TRACKDB_HH_

//----------------------------------------------------------------------------//
//              end of file TrackDB.hh
//----------------------------------------------------------------------------//
