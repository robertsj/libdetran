//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   TrackDB.hh
 * \brief  TrackDB 
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef TRACKDB_HH_
#define TRACKDB_HH_

#include "geometry/geometry_export.hh"
#include "Track.hh"
#include "angle/QuadratureMOC.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_geometry
{

/*!
 *  \class TrackDB
 *  \brief Database of tracks.
 *
 *  The track database contains all the information needed
 *  to represent the MOC problem with respect to space
 *  and angle.
 *
 *  Currently, this is limited to 2-D tracking.  Tracks
 *  are only assigned for azimuths over [0, pi], using
 *  symmetry for [pi, 2*pi].
 *
 *  Tracks are stored in order in order of decreasing
 *  y and increasing x for 0, pi/2 and decreasing
 *  x for pi/2, pi.  This is how tracks are indexed.
 *  The track index in the other two octancts keeps
 *  the index of their reflection.
 *
 */
/*!
 *  \example geometry/test/test_TrackDB.cc
 *
 *  Test of TrackDB class
 */
class GEOMETRY_EXPORT TrackDB
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<TrackDB>               SP_trackdb;
  typedef detran_angle::QuadratureMOC::SP_quadrature  SP_quadrature;
  typedef Track::SP_track                             SP_track;
  typedef std::vector<SP_track>                       vec_track;
  typedef std::vector<vec_track>                      vec2_track;
  typedef detran_utilities::size_t                    size_t;
  typedef detran_utilities::vec_int                   vec_int;
  typedef detran_utilities::vec_dbl                   vec_dbl;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param num_azimuths   Number of azimuths in first two octants
   */
  TrackDB(int num_azimuths, int num_regions, SP_quadrature quad)
    : d_number_azimuths(num_azimuths)
    , d_number_regions(num_regions)
    , d_quadrature(quad)
    , d_tracks(num_azimuths)
    , d_cos_phi(num_azimuths, 0.0)
    , d_sin_phi(num_azimuths, 0.0)
    , d_spacing(num_azimuths, 0.0)
  {
    Require(d_number_azimuths > 0);
    Require(d_number_regions > 0);
  }

  ~TrackDB(){}

  /// \name Getters
  /// \{

  SP_track track(size_t a, size_t t)
  {
    Require(a < d_tracks.size());
    Require(t < d_tracks[a].size());
    return d_tracks[a][t];
  }

  int number_tracks_angle(size_t a) const
  {
    Require(a < d_tracks.size());
    return d_tracks[a].size();
  }

  int number_angles() const
  {
    return d_tracks.size();
  }

  double spacing(size_t a) const
  {
    Require(a < d_spacing.size());
    return d_spacing[a];
  }

  /// \}

  void add_track(size_t a, SP_track t)
  {
    Require(a < d_tracks.size());
    d_tracks[a].push_back(t);
  }

  void setup_angle(size_t a, double c_phi, double s_phi, double space)
  {
    Require(a < d_tracks.size());
    Require(space > 0.0);
    d_cos_phi[a] = c_phi;
    d_sin_phi[a] = s_phi;
    d_spacing[a] = space;
  }

  /// Normalize the tracks given a vector of true volumes.
  void normalize(vec_dbl &volume);

  /// Pretty display of all track
  void display() const;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Quadrature
  SP_quadrature d_quadrature;
  /// Number of azimuths in first two octants
  int d_number_azimuths;
  /// Number of flat source regions
  int d_number_regions;
  /// Tracks by [azimuth][space]
  vec2_track d_tracks;
  /// Azimuthal cosines.
  vec_dbl d_cos_phi;
  /// Azimuthal sines.
  vec_dbl d_sin_phi;
  /// Constant track spacing for each angle.
  vec_dbl d_spacing;

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<TrackDB>)

} // end namespace detran_geometry

#endif // TRACKDB_HH_ 

//---------------------------------------------------------------------------//
//              end of file TrackDB.hh
//---------------------------------------------------------------------------//
