//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   TransportManager.hh
 *  @brief  TransportManager
 *  @author Jeremy Roberts
 *  @date   Nov 6, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_TRANSPORTMANAGER_HH_
#define detran_TRANSPORTMANAGER_HH_

#include "callow/utils/Initialization.hh"
#include "utilities/InputDB.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"

namespace detran
{

/**
 *  @class TransportManager
 *  @brief Base class for all transport managers.
 *
 *  The purpose of an untemplated base mesh is to simplify
 *  initialization of libraries, etc.
 */
class TransportManager
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATION
  //-------------------------------------------------------------------------//

  // Basic spatial discretization categories
  enum EQTYPES
  {
    SN, MOC, DIFF, END_EQTYPES
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input           SP_input;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_material::Material::SP_material        SP_material;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor.  This initializes all external libraries.
  TransportManager(int argc, char** argv)
    : d_flag(true)
  {
    callow_initialize(argc, argv);
  }

  /// Default constructor
  TransportManager()
    : d_flag(false)
  {
    /* ... */
  }

  /// Destructor.  This initializes all external libraries.
  virtual ~TransportManager()
  {
    if (d_flag) callow_finalize();
  }

private:

  /// Flag for initialization.
  bool d_flag;
  /// User input
  SP_input d_input;
  /// Material database
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;

};


} // end namespace detran

#endif // TRANSPORTMANAGER_HH_ 

//---------------------------------------------------------------------------//
//              end of file TransportManager.hh
//---------------------------------------------------------------------------//
