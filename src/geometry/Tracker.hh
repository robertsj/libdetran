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
#include "MeshMOC.hh"
#include "angle/QuadratureMOC.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_geometry
{

/**
 *  @class Tracker
 *  @brief Track a mesh
 */
class GEOMETRY_EXPORT Tracker
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Tracker>               SP_tracker;
  typedef Mesh::SP_mesh                               SP_mesh;
  typedef detran_angle::QuadratureMOC::SP_quadrature  SP_quadrature;
  typedef Track::SP_track                             SP_track;
  typedef TrackDB::SP_trackdb                         SP_trackdb;
  typedef detran_utilities::vec_dbl                   vec_dbl;
  typedef detran_utilities::Point                     Point;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    mesh            mesh to track
   *  @param    quadrature      product quadrature
   */
  Tracker(SP_mesh mesh, SP_quadrature quadrature);

  /// SP constructor
  static SP_tracker Create(SP_mesh mesh, SP_quadrature quadrature);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  SP_trackdb trackdb() const
  {
    Require(d_trackdb);
    return d_trackdb;
  }

  SP_mesh meshmoc()
  {
    // Create the MOC mesh
    SP_mesh newmesh(new MeshMOC(d_mesh, d_trackdb));
    return newmesh;
  }

  // Normalize the track segments based on actual volumes.
  void normalize();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_mesh d_mesh;
  SP_quadrature d_quadrature;
  SP_trackdb d_trackdb;
  // Number azimuths per octant
  int d_number_azimuths;
  vec_dbl d_x;
  vec_dbl d_y;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  void generate_tracks();
  void find_starting_cell(Point enter, double tan_phi, int *IJ);

};

} // end namespace detran_geometry

#endif /* detran_geometry_TRACKER_HH_ */

//----------------------------------------------------------------------------//
//              end of file Tracker.hh
//----------------------------------------------------------------------------//
