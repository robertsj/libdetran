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
//#include "detran_config.h"
//
//#ifdef DETRAN_ENABLE_SILO

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"

// Utilities
#include "InputDB.hh"

// System
#include "silo.h"

namespace detran_postprocess
{

/*!
 *  \class SiloOutput
 *  \brief Write flux data to a silo file.
 */
class SiloOutput
{

public:

  typedef detran::InputDB::SP_input     SP_input;
  typedef detran::Material::SP_material SP_material;
  typedef detran::Mesh::SP_mesh         SP_mesh;
  typedef detran::State::SP_state       SP_state;

  /*!
   *  \brief Constructor
   *  \param input    Input database
   *  \param mesh     Cartesian mesh
   */
  SiloOutput(SP_input input, SP_mesh mesh);

  /// Destructor
  ~SiloOutput();

  /// Initialize file
  void initialize();

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

  /*!
   *  \brief Write the multigroup scalar flux to file.
   *  \param state  State vector container
   */
  void write_flux(SP_state state);

private:

  /// Input database
  SP_input d_input;

  /// Problem mesh
  SP_mesh d_mesh;

  /// Material database
  SP_material d_material;

  /// Silo file.
  DBfile *d_silofile;

  /// The Silo file is created
  bool d_initialized;

};

} // end namespace detran_postprocess

//#endif // DETRAN_ENABLE_SILO

#endif // SILOOUTPUT_HH_ 

//---------------------------------------------------------------------------//
//              end of file SiloOutput.hh
//---------------------------------------------------------------------------//
