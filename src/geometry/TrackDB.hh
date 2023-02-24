//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TrackDB.hh
 *  @brief TrackDB class definition
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
#include <memory>
#include <vector>

namespace detran_geometry
{

/**
 *  @class TrackDB
 *  @brief Database of tracks.
 *
 *  NEW IDEA
 *   - Have tracker add all tracks to a single vector of tracks.
 *   - Construct TrackDB from this vector
 *   - Provide groupings that make sense (e.g., tracks for an octant, tracks for a face, etc.)
 *
 *  Tracks are defined according to an azimuth, polar angle (3-D), and
 *  track spacing by Tracker.  Within TrackDB, a vector of Track
 *  pointers is stored as an element of a 2-D vector indexed by
 *  azimuth and polar angle.
 *
 *  In 2-D, tracks are stored only for the first two quadrants.
 *  In 3-D, tracks are similarly stored only for the first four octants.
 *
 *  For the other quadrants and octants, track reverse iterat
 *
 *  Tracks are not guaranteed to be stored in any particular order within
 *  an angle.  The client needs to post process to obtain proper indexing,
 *  e.g. for boundary conditions.
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

  typedef std::shared_ptr<detran_angle::ProductQuadrature>    SP_quadrature;
  typedef std::vector<std::shared_ptr<Track>>                 vec_track;
  typedef std::map<std::pair<int, int>, vec_track>            map_track;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param q        product quadrature
   *  @param tracks   vector of track pointers
   *  @param symmetry is symmetry used to reduce the number of tracks by half?
   */
  TrackDB(SP_quadrature q);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Get the tracks.
  const map_track& tracks() const
  {
    return d_tracks;
  };
  map_track& tracks()
  {
    return d_tracks;
  };

  /// Normalize the tracks given a vector of true region volumes.
  void normalize(const std::vector<double> &volume);
  /// Sort the tracks from left-to-right, bottom-to-top on an incident side.
  void sort();
  /// Pretty display of all tracks.
  void display() const;

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//
  /// map of octant,angle indices to vector of tracks.
  map_track d_tracks;
  /// Quadrature
  SP_quadrature d_quadrature;


};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<TrackDB>)

} // end namespace detran_geometry

#endif // detran_geometry_TRACKDB_HH_

//----------------------------------------------------------------------------//
//              end of file TrackDB.hh
//----------------------------------------------------------------------------//
