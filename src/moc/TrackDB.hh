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

// Detran
#include "Track.hh"
#include "QuadratureMOC.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

#include <vector>

namespace detran
{

/*!
 *  \class TrackDB
 *  \brief Database of tracks.
 *
 *  The track database contains all the information needed
 *  to represent the MOC problem with respect to space
 *  and angle.  Essentially, it contains the following:
 *
 */
class TrackDB: public Object
{

public:

  typedef SP<TrackDB>                   SP_trackdb;
  typedef QuadratureMOC::SP_quadrature  SP_quadrature;
  typedef Track::SP_track               SP_track;
  typedef std::vector<SP_track>         vec_track;
  typedef std::vector<vec_track>        vec2_track;

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

  SP_track track(u_int a, u_int t)
  {
    Require(a < d_tracks.size());
    Require(t < d_tracks[a].size());
    return d_tracks[a][t];
  }

  int number_tracks_angle(u_int a) const
  {
    Require(a < d_tracks.size());
    return d_tracks[a].size();
  }

  int number_angles() const
  {
    return d_tracks.size();
  }

  double spacing(u_int a) const
  {
    Require(a < d_spacing.size());
    return d_spacing[a];
  }

  /// \}

  void add_track(u_int a, SP_track t)
  {
    Require(a < d_tracks.size());
    d_tracks[a].push_back(t);
  }

  void setup_angle(u_int a, double c_phi, double s_phi, double space)
  {
    Require(a < d_tracks.size());
    Require(space > 0.0);
    d_cos_phi[a] = c_phi;
    d_sin_phi[a] = s_phi;
    d_spacing[a] = space;
  }

  /// Normalize the tracks given a vector of true volumes.
  void normalize(vec_dbl &volume);

  void display() const;

  bool is_valid() const
  {
    return true;
  }

private:

  /// \name Private Data
  /// \{

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

  /// \}


};


} // end namespace detran

#endif // TRACKDB_HH_ 

//---------------------------------------------------------------------------//
//              end of file TrackDB.hh
//---------------------------------------------------------------------------//
