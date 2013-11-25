//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SiloOutput.hh
 *  @brief  SiloOutput
 *  @author Jeremy Roberts
 *  @date   Jul 27, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_ioutils_SILOOUTPUT_HH_
#define detran_ioutils_SILOOUTPUT_HH_

#include "ioutils/ioutils_export.hh"
#include "detran_config.hh"
#include "angle/Quadrature.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "transport/State.hh"
#include "utilities/InputDB.hh"
#ifdef DETRAN_ENABLE_SILO
#include "silo.h"
#endif

namespace detran_ioutils
{

/**
 *  @class SiloOutput
 *  @brief Write data to a Silo file.
 *
 *  This class wraps the Silo API to allow clients to write mesh data to
 *  binary.  Special functions automate writing multigroup scalar and
 *  angular fluxes.  Additionally, time-dependent fluxes can be written.
 */
/**
 *  @example ioutils/test_SiloOutput.cc
 *  @brief   Test of SiloOutput
 */
class IOUTILS_EXPORT SiloOutput
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<SiloOutput>        SP_silooutput;
  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_material::Material::SP_material  SP_material;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran::State::SP_state                 SP_state;
  typedef detran_angle::Quadrature::SP_quadrature SP_quadrature;
  typedef detran_utilities::vec_dbl               vec_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input    Input database
   *  @param mesh     Cartesian mesh
   */
  SiloOutput(SP_mesh mesh);

  /// Destructor
  ~SiloOutput();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Initialize file
  bool initialize(const std::string filename);

  /**
   *  @brief Write a mesh map to file.
   *  @param key  Key of the map to write
   *  @return     True for successful write
   */
  bool write_mesh_map(const std::string &key);

  /**
   *  @brief Write the multigroup scalar flux moments to file.
   *  @param state  State vector container
   *  @return     True for successful write
   */
  bool write_scalar_flux(SP_state state);

  /**
   *  @brief Write the angular flux to file.
   *  @param state  State vector container
   *  @return       True for successful write
   */
  bool write_angular_flux(SP_state state, SP_quadrature quad);

  /**
   *  @brief Write the time dependent fluxes to file
   *  @param  step   Time step
   *  @param  state  State vector container
   *  @return        True for successful write
   */
  bool write_time_flux(const int step, SP_state state, bool do_psi);

  /// Write a mesh scalar field
  bool write_scalar_field(const std::string &key, const vec_dbl &data);

  /// Write a mesh vector field
  bool write_vector_field(const std::string &key,
                          const vec_dbl     &data_i,
                          const vec_dbl     &data_j = vec_dbl(0),
                          const vec_dbl     &data_k = vec_dbl(0));

  /// Close file
  void finalize();

  /**
   *  @brief Add the material database
   *
   *  The user can create maps of the cross section
   *  values.
   *
   *  @param material Material database
   */
  void add_material(SP_material material)
  {
    Require(material);
    d_material = material;
  }

  bool make_directory(const std::string &dir);

  bool set_directory(const std::string &dir);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Problem mesh
  SP_mesh d_mesh;
  /// Material database
  SP_material d_material;
  /// Silo file.
#ifdef DETRAN_ENABLE_SILO
  DBfile *d_silofile;
#endif
  /// The Silo file is created
  bool d_initialized;
  /// Problem dimension
  int d_dimension;
  /// Direction dimensions
  int d_dims[3];

};

IOUTILS_TEMPLATE_EXPORT(detran_utilities::SP<SiloOutput>)

} // end namespace detran_ioutils

#endif // detran_ioutils_SILOOUTPUT_HH_

//---------------------------------------------------------------------------//
//              end of file SiloOutput.hh
//---------------------------------------------------------------------------//
