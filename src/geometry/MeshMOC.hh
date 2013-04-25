//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MeshMOC.hh
 * \brief  MeshMOC 
 * \author Jeremy Roberts
 * \date   Jun 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MESHMOC_HH_
#define MESHMOC_HH_

#include "geometry/geometry_export.hh"
#include "Mesh.hh"
#include "TrackDB.hh"

namespace detran_geometry
{

/*!
 *  \class MeshMOC
 *  \brief Mesh with tracking information.
 *  \todo  This is bloat.  There should be one mesh object, or
 *         perhaps one geometry object that contains a mesh and tracking
 */
class GEOMETRY_EXPORT MeshMOC: public Mesh
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<MeshMOC>     SP_mesh;
  typedef Mesh                              Base;
  typedef Base::SP_mesh                     SP_base;
  typedef TrackDB::SP_trackdb               SP_trackdb;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  MeshMOC(SP_base mesh, SP_trackdb tracks)
    : Mesh(*mesh)
    , d_tracks(tracks)
  {
    Require(d_tracks);
  }

  SP_trackdb tracks() const
  {
    return d_tracks;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Track database
  SP_trackdb d_tracks;

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<MeshMOC>)

} // end namespace detran_geometry

#endif // MESHMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file MeshMOC.hh
//---------------------------------------------------------------------------//
