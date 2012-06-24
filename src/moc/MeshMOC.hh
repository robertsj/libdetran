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

// Detran
#include "Mesh.hh"
#include "TrackDB.hh"


namespace detran
{

class MeshMOC: public Mesh
{

public:

  typedef SP<MeshMOC>           SP_mesh;
  typedef Mesh                  Base;
  typedef Base::SP_mesh         SP_base;
  typedef TrackDB::SP_trackdb   SP_trackdb;

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

  /// \name Private Data
  /// \{

  /// Track database
  SP_trackdb d_tracks;

  /// \}

};

} // end namespace detran

#endif // MESHMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file MeshMOC.hh
//---------------------------------------------------------------------------//
