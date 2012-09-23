//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SiloOutput.hh
 * \brief  SiloOutput 
 * \author Jeremy Roberts
 * \date   Jul 27, 2012
 */
//---------------------------------------------------------------------------//

#ifndef SILOOUTPUT_HH_
#define SILOOUTPUT_HH_

#include "detran_config.hh"
#include "angle/Quadrature.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/State.hh"
#include "utilities/InputDB.hh"
#include "silo.h"

namespace detran_ioutils
{

/*!
 *  \class SiloOutput
 *  \brief Write mesh data to a Silo file.
 */
class SiloOutput
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_material::Material::SP_material  SP_material;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran::State::SP_state                 SP_state;
  typedef detran_angle::Quadrature::SP_quadrature SP_quadrature;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param input    Input database
   *  \param mesh     Cartesian mesh
   */
  SiloOutput(SP_mesh mesh);

  /// Destructor
  ~SiloOutput();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Initialize file
  bool initialize(const std::string filename);

  /*!
   *  \brief Write a mesh map to file.
   *  \param key  Key of the map to write
   *  \return     True for successful write
   */
  bool write_mesh_map(const std::string &key);

  /*!
   *  \brief Write the multigroup scalar flux moments to file.
   *  \param state  State vector container
   *  \return     True for successful write
   */
  bool write_scalar_flux(SP_state state);

  /*!
   *  \brief Write the angular flux to file.
   *  \param state  State vector container
   *  \return       True for successful write
   */
  bool write_angular_flux(SP_state state, SP_quadrature quad);

  /// Close file
  void finalize();

  /*!
   *  \brief Add the material database
   *
   *  The user can create maps of the cross section
   *  values.
   *
   *  \param material Material database
   */
  void add_material(SP_material material)
  {
    Require(material);
    d_material = material;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Problem mesh
  SP_mesh d_mesh;

  /// Material database
  SP_material d_material;

  /// Silo file.
  DBfile *d_silofile;

  /// The Silo file is created
  bool d_initialized;

  /// Problem dimension
  int d_dimension;

  /// Direction dimensions
  int d_dims[3];

};

} // end namespace detran_ioutils

#endif // SILOOUTPUT_HH_ 

//---------------------------------------------------------------------------//
//              end of file SiloOutput.hh
//---------------------------------------------------------------------------//
