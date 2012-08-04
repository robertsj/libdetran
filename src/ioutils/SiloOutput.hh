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

// Configuration
#include "detran_config.h"

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"

// Utilities
#include "InputDB.hh"

// System
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

  typedef detran::InputDB::SP_input         SP_input;
  typedef detran::Material::SP_material     SP_material;
  typedef detran::Mesh::SP_mesh             SP_mesh;
  typedef detran::State::SP_state           SP_state;
  typedef detran::Quadrature::SP_quadrature SP_quadrature;

  /*!
   *  \brief Constructor
   *  \param input    Input database
   *  \param mesh     Cartesian mesh
   */
  SiloOutput(SP_mesh mesh);

  /// Destructor
  ~SiloOutput();

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
