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
 */
class TrackDB: public Object
{

public:

  typedef SP<TrackDB>             SP_trackdb;
  typedef Track::SP_track         SP_track;
  typedef std::vector<SP_track>   vec_track;
  typedef std::vector<vec_track>  vec2_track;

  TrackDB(int num_azimuth)
    : d_tracks(num_azimuth)
  { /* ... */ }

  ~TrackDB(){}

  SP_track track(u_int a, u_int t)
  {
    Require(a < d_tracks.size());
    Require(t < d_tracks[a].size());
    return d_tracks[a][t];
  }

  void add_track(u_int a, SP_track t)
  {
    Require(a < d_tracks.size());
    d_tracks[a].push_back(t);
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

  void display() const;

  bool is_valid() const
  {
    return true;
  }

private:

  /// \name Private Data
  /// \{

  /// Tracks by [azimuth][space]
  vec2_track d_tracks;

  /// \}


};


} // end namespace detran

#endif // TRACKDB_HH_ 

//---------------------------------------------------------------------------//
//              end of file TrackDB.hh
//---------------------------------------------------------------------------//
