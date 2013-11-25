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
#include "utilities/Iterators.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_geometry
{

/**
 *  @class TrackDB
 *  @brief Database of tracks.
 *
 *  The track database contains all the information needed to represent the
 *  MOC problem with respect to space and angle.  In 2-D, tracks are
 *  stored only for the first two quadrants, phi = [0, pi].  In 3-D, tracks are
 *  similarly stored only for the first four octants, phi = [0, 2*pi] and
 *  theta = [0, pi/2].
 *
 *  Tracks are not guaranteed to be stored in any particular order within
 *  an angle.  The client needs to post process to obtain proper indexing,
 *  e.g. for boundary conditions.
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
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum DIRECTIONS
  {
    BACKWARD, FORWARD, END_DIRECTIONS
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<TrackDB>                       SP_trackdb;
  typedef detran_angle::ProductQuadrature::SP_quadrature      SP_quadrature;
  typedef Track::SP_track                                     SP_track;
  typedef std::vector<SP_track>                               vec_track;
  typedef std::vector<vec_track>                              vec2_track;
  typedef std::vector<vec2_track>                             vec3_track;
  typedef detran_utilities::size_t                            size_t;
  typedef const size_t                                        c_size_t;
  typedef detran_utilities::vec_int                           vec_int;
  typedef detran_utilities::vec_dbl                           vec_dbl;

  typedef vec_track::iterator                                 iterator_angle;
  //typedef detran_utilities::Reversible<vec_track>         iterator;


  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    q   product quadrature
   */
  TrackDB(SP_quadrature q);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Get a track
   *  @param    a       azimuth index
   *  @param    p       polar index
   *  @param    t       track index within [a, p]
   */
  SP_track track(c_size_t a, c_size_t p, c_size_t t);

  struct iterator_anglea
  {

  };


  //@{
  /// Iterators to all tracks
//  iterator begin(bool forward = true);
//  iterator end(bool forward = true);
  /// Iterators to tracks for a given angle
  iterator_angle begin(c_size_t a, c_size_t p);
  iterator_angle end(c_size_t a, c_size_t p);
  //@}

  /**
   *  @brief Number of tracks in a given angle
   *  @param    a       azimuth index
   *  @param    p       polar index
   */
  size_t number_tracks(c_size_t a, c_size_t p = 0) const;

  /**
   *  @brief Add a track
   *  @param    a       azimuth index
   *  @param    p       polar index
   *  @param    t       track index within azimuths
   */
  void add_track(c_size_t a, c_size_t p, SP_track t);

  /// Normalize the tracks given a vector of true region volumes.
  void normalize(const vec_dbl &volume);

  /// Get the number of azimuths stored
  size_t number_azimuths() const;

  /// Get the number of polar angles stored
  size_t number_polar() const;

  /// Return the dimension
  size_t dimension() const;

  /// Sort the tracks from left-to-right, bottom-to-top on an incident side
  void sort();

  /// Pretty display of all track
  void display() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Quadrature
  SP_quadrature d_quadrature;
  /// Dimension
  size_t d_dimension;
  /// Number of azimuths stored
  size_t d_number_azimuths;
  /// Number of polar angles stored
  size_t d_number_polar;
  /// Tracks by [azimuth][polar][space]
  vec3_track d_tracks;

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<TrackDB>)

} // end namespace detran_geometry

#endif // detran_geometry_TRACKDB_HH_

//----------------------------------------------------------------------------//
//              end of file TrackDB.hh
//----------------------------------------------------------------------------//
